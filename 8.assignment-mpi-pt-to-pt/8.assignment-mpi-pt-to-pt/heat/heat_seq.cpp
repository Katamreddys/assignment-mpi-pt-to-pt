#include <iostream>
#include <chrono>
#include <math.h>


#define DEBUG 0 

#ifdef __cplusplus
extern "C" {
#endif
  double generate2DHeat(long n, long global_i, long global_j);
  int check2DHeat(long n, long global_i, long global_j, double v, long k); //this function return 1 on correct. 0 on incorrect. Note that it may return 1 on incorrect. But a return of 0 means it is definitely incorrect

#ifdef __cplusplus
}
#endif



/**
 *
 */
void showHeat(double **H, int m, int n) {
  for (int i=0; i<m; ++i) {
    for (int j=0; j<n; ++j) {
      std::cout<<H[i][j]<<" ";
    }
    std::cout<<std::endl;
  }
}


/**
 *
 */
double** prev_heat(double ***sub_heat, int k) {
  return sub_heat[(k-1)%2];
}


/**
 *
 */
double** curr_heat(double ***sub_heat, int k) {
  return sub_heat[k%2];
}


/**
 * block[9] 
 *
 * 0 1 2
 * 3 4 5   ~~~ (i, j) is position 4
 * 6 7 8
 *
 */
void heat_eq(double **prev, double **curr, int i, int j, int N) {

  int stencil = 5;
  double *block = new double[stencil];
  // init all to center of block, prev[i][j] 
  for (int t=0; t<stencil; ++t) {
    block[t] = prev[i][j]; 
  }
  
  // column LEFT of center
  if ((j-1) >= 0) {
    block[1] = prev[i][j-1];
  }
  
  // above center
  if ((i-1) >= 0) 
    block[0] = prev[i-1][j];
  // below center
  if ((i+1) < N) 
    block[4] = prev[i+1][j];
  
  // column RIGHT of center
  if ((j+1) < N) {
    block[3] = prev[i][j+1];
  }
  
#if DEBUG //
  std::cout<<"i: "<<i;
  std::cout<<", j: "<<j<<" = [ ";
  for (int t=0; t<9; ++t) {
    if (t % 3 == 0) std::cout<<std::endl;
    std::cout<<block[t]<<" ";
  }
  std::cout<<std::endl<<"]"<<std::endl<<std::endl;
#endif

  // reset center of block i,j
  curr[i][j] = 0.0f; 
  // sum up block
  for (int t=0; t<stencil; ++t) {
    curr[i][j] += block[t];
  }
  // average
  curr[i][j] /= (float) stencil;

  delete[] block;
}


/**
 *
 */
int main(int argc, char* argv[]) {

  if (argc < 5) {
    std::cerr<<"usage: "<<argv[0]<<" <N> <K> <rank> <P>"<<std::endl;
    return -1;
  }

  // declare and init command line params
  long N, K, rank, P;
  N = std::atol(argv[1]);
  K = std::atol(argv[2]);
  rank = std::atol(argv[3]);
  P = std::atol(argv[4]);


  double ***sub_heat = new double**[2];

  int sub_N = N/sqrt(P);

  for (int k=0; k<2; ++k) {
    sub_heat[k] = new double*[sub_N];
    for (int i=0; i<sub_N; ++i) { // dont init buffers, i==0, i==m-1
      sub_heat[k][i] = new double[sub_N];
    }
#if DEBUG
    if (k == 0) showHeat(sub_heat[k], N, N);
#endif
  }

  long row = rank%lround(sqrt(P));
  long col = rank/lround(sqrt(P));

  std::cout<<"row: "<<row<<" col: "<<col<<std::endl;
  
  // generate 2D array
  //generate2DHeat(sub_heat[0], N, rank, P);
  {
    double** h = sub_heat[0];
    for (long i = 0; i<sub_N; ++i)
      for (long j = 0; j<sub_N; ++j) {
	long global_i = i + sub_N * row;
	long global_j = j + sub_N * col;
	h[i][j] = generate2DHeat(N, global_i, global_j);
      }
	  
  }

  showHeat(sub_heat[0], sub_N, sub_N);

  // convient handles for Heat^k-1, Heat^k
  double **prev, **curr;

  // begin timing
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();


  // repeat K times
  for (int k=1; k<=K; ++k) {

    prev = prev_heat(sub_heat, k); // Heat^k-1
    curr = curr_heat(sub_heat, k); // Heat^k

    // compute heat equation using owned data
    for (int i=0; i<sub_N; ++i) { 
      for (int j=0; j<sub_N; ++j) {
        heat_eq(prev, curr, i, j, sub_N);
      }
    }
  
   {
     double** h = curr;
     bool correct = true;
     for (long i = 0; i<sub_N; ++i)
       for (long j = 0; j<sub_N; ++j) {
	 long global_i = i + sub_N * row;
	 long global_j = j + sub_N * col;
	 correct &= check2DHeat(N, global_i, global_j, h[i][j], k);
	 if (!correct) {
	   std::cerr<<"incorrect computation at iteration ("<<global_i<<","<<global_j<<") "<<k<<" new check"<<std::endl;
	 }

       }

  }
  
    std::cout<<"_________________________"<<std::endl;
    showHeat(curr, sub_N, sub_N);

#if DEBUG
    showHeat(curr, sub_N, sub_N);
    std::cout<<"_________________________"<<std::endl;
#endif
  }

  // end timing
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  // report time
  std::cerr<<elapsed_seconds.count()<<std::endl;

  return 0;
}

