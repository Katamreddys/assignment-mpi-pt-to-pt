#include <mpi.h>
#include <math.h>
#include <iostream>

using namespace std;

#ifdef __cplusplus
extern "C" {
 #endif

  int check2DHeat(double** H, long n, long rank, long P, long k); //this assumes array of array and grid block decomposition

 #ifdef __cplusplus
}
#endif

/***********************************************
 *         NOTES on check2DHeat.
 ***********************************************
 *         
 *  First of, I apologize its wonky. 
 *
 *  Email me ktibbett@uncc.edu with any issues/concerns with this. Dr. Saule or the other
 *    TA's are not familiar with how it works. 
 *
 * Params:
 *  n - is the same N from the command line, NOT the process's part of N
 *  P - the total amount of processes ie what MPI_Comm_size gives you.
 *  k - assumes n/2 > k-1 , otherwise may return false negatives.
 *
 *   
 * Disclaimer:
 ***
 *** Broken for P is 9. Gives false negatives, for me it was always
 ***  ranks 0, 3, 6. I have not found issues with 1, 4, or 16, and these
 ***  are what `make test` will use.
 ***
 *
 * Usage:
 *  When code is WRONG returns TRUE. Short example below
 *  if (check2DHeat(...)) {
 *    // oh no it is (maybe) wrong  
 *    std::cout<<"rank: "<<rank<<" is incorrect"<<std::endl;
 *  }
 *
 *
 *
 *  I suggest commenting this out when running the bench
 *
 *
 * - Kyle
 *
 *************/


// Use similarily as the genA, genx from matmult assignment.
double genH0(long row, long col, long n) {
  double val = (double)(col == (n/2));
  return val;
}



int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1;
  }

  // declare and init command line params
  MPI_Init(&argc,&argv);
  long n, K;
  n = atol(argv[1]);
  K = atol(argv[2]);

  int worldrank,np;
   
   MPI_Comm_rank(MPI_COMM_WORLD,&worldrank);
   MPI_Comm_size(MPI_COMM_WORLD,&np);

   int p = sqrt(np);
   long each_div = n/p; 
   
   int row_div = worldrank/p,col_div = worldrank%p;
  // use double for heat 2d 

   double**  H = new double*[each_div];
   double**  G = new double*[each_div];
   for(long i=0;i<each_div;i++)
     {
     H[i] = new double[each_div];
     G[i] = new double[each_div];
     
   }
   
   long rowstart = (row_div*each_div),colstart = (col_div*each_div);
   long rowend = rowstart+each_div,colend = colstart+each_div;

  for (long row = rowstart,rset=0; row<rowend; row++,rset++) {
    for (long col= colstart,cset=0; col<colend; col++,cset++) {
       H[rset][cset] = genH0(row, col,n);
    }
  } 
  // write code here
  
  double *tp = new double[each_div];
  double *btm = new double[each_div];
  double *rgt = new double[each_div];
  double *lft = new double[each_div];
  double *tp1 = new double[each_div];
  double *btm1 = new double[each_div];
  double *rgt1 = new double[each_div];
  double *lft1 = new double[each_div];
  
  int top,bottom,left,right;
  left =col_div?worldrank-1:-1;
  right = (col_div == (p-1))?-1:worldrank+1;
  top = worldrank-p;
  bottom = worldrank+p;
  MPI_Status status[4];
  MPI_Request requests[8];
  long count = 0;
  double start = MPI_Wtime();
  for (long it = 0; it<K; it++) 
  {
     for(long ind =0,ite = 0;ind < each_div;ind++)
       {
   lft[ite] = H[ind][0];
   rgt[ite] = H[ind][each_div-1];
   tp[ite]  = H[0][ind];
   btm[ite] = H[each_div-1][ind];
   lft1[ite] = H[ind][0];
   rgt1[ite] = H[ind][each_div-1];
   tp1[ite]  = H[0][ind];
   btm1[ite] = H[each_div-1][ind];
   ite++;
       }
     count = 0;
     //  cout<<"I am rank "<<worldrank<<endl; 
     if(top>=0){
       MPI_Isend(tp1,each_div,MPI_DOUBLE,top,0,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(right != -1){
       MPI_Isend(rgt1,each_div,MPI_DOUBLE,right,1,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(bottom < np){
       MPI_Isend(btm1,each_div,MPI_DOUBLE,bottom,2,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     if(left != -1){
       MPI_Isend(lft1,each_div,MPI_DOUBLE,left,3,MPI_COMM_WORLD,&requests[count]);
       count++;
     }
     for(long i=1;(i+1)<each_div;i++)
       {
   for(long j=1;(j+1)<each_div;j++)
     {
       G[i][j] = (H[i][j] + H[i-1][j] + H[i+1][j] + H[i][j-1] + H[i][j+1])/5;
     }
       }
     // cout<<"I am rank "<<worldrank<<" and sent "<<count<<" messages"<<endl; 
     // MPI_Waitall(count,requests,status);
     if(top >= 0){
       MPI_Recv(tp,each_div,MPI_DOUBLE,top,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(right != -1){
       MPI_Recv(rgt,each_div,MPI_DOUBLE,right,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      
     }
     if(bottom < np){
       MPI_Recv(btm,each_div,MPI_DOUBLE,bottom,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     if(left != -1){
       MPI_Recv(lft,each_div,MPI_DOUBLE,left,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       
     }
     // cout<<"I am rank "<<worldrank<<" and received "<<count<<" messages"<<endl;
     G[0][0] = (H[0][0] + tp[0] + lft[0] + H[0][1] + H[1][0])/5;
     G[0][each_div-1] = (H[0][each_div-1] + tp[each_div-1] + rgt[0]+H[0][each_div-2]+H[1][each_div-1])/5;
     G[each_div-1][0] = (H[each_div-1][0] + btm[0] + lft[each_div-1]+H[each_div-1][1]+H[each_div-2][0])/5;
     G[each_div-1][each_div-1] = (H[each_div-1][each_div-1] + btm[each_div-1] +
          rgt[each_div-1]+H[each_div-1][each_div-2]+H[each_div-2][each_div-1])/5;
     for(long i=1,j=1;i<each_div-1;i++,j++)
   {
      G[0][i] = (H[0][i]+H[1][i]+tp[i]+H[0][i-1]+H[0][i+1])/5;
      G[each_div-1][i] = (H[each_div-1][i]+H[each_div-2][i]+btm[i]+H[each_div-1][i-1]+H[each_div-1][i+1])/5;
      G[j][0] = (H[j][0]+lft[j]+H[j][1]+H[j-1][0]+H[j+1][0])/5;
      G[j][each_div-1] = (H[j][each_div-1]+H[j-1][each_div-1]+H[j+1][each_div]+rgt[j]+H[j][each_div-2])/5;
   }
      H = G;
      // check2DHeat(double** H, long n, long rank, long P, long k)
      check2DHeat(H,n,worldrank,np,it);
      MPI_Waitall(count,requests,status);
             
   }
  if(worldrank == 0)
    {
       double end = MPI_Wtime();
       cerr<<end-start<<endl; 
    }
 for(long i=0;i<each_div;i++)
   delete[] H[i];
  delete[] H;
  delete[] tp;
  delete[] btm;
  delete[] rgt;
  delete[] lft;
  MPI_Finalize();

  return 0;
}