#include <mpi.h>
#include <math.h>
#include <iostream>
#include <chrono>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

  double generate2DHeat(long n, long global_i, long global_j);

  int check2DHeat(double** H, long n, long rank, long P, long k); 

#ifdef __cplusplus
}
#endif

void calculate_2d_heat(long block, double** Current, double** Previous, double* move_left, double* move_right, double* move_up, double* move_down){

  int row, col;
  
  for (row = 0; row < block; ++row) { 
    for (col = 0; col < block; ++col) {
      if(row == 0){
    Current[row][col] = (move_up[col] + Previous[row][col-1] + Previous[row][col] + Previous[row][col+1] + Previous[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1){
    Current[row][col] = (move_down[col] + Previous[row-1][col] + Previous[row][col-1] + Previous[row][col] + Previous[row][col+1])/(static_cast<double>(5));
      }
      else if(col == 0){
    Current[row][col] = (move_left[row] + Previous[row-1][col] + Previous[row][col] + Previous[row][col+1] + Previous[row+1][col])/(static_cast<double>(5));
      }
      else if(col == block-1){
    Current[row][col] = (move_right[row] + Previous[row-1][col] + Previous[row][col-1] + Previous[row][col]+ Previous[row+1][col])/(static_cast<double>(5));
      }
      else if(row == 0 && col == 0){
    Current[row][col] = (move_left[row] + move_up[col] + Previous[row][col] + Previous[row][col+1] + Previous[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1 && col == 0){
    Current[row][col] = (move_left[row] + move_down[col] + Previous[row-1][col] + Previous[row][col] + Previous[row][col+1])/(static_cast<double>(5));
      }
      else if(row == 0 && col == block-1){
    Current[row][col] = (move_right[row] + move_up[col] + Previous[row][col-1] + Previous[row][col] + Previous[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1 && col == block-1){
    Current[row][col] = (move_right[row] + move_down[col] + Previous[row-1][col] + Previous[row][col-1] + Previous[row][col])/(static_cast<double>(5));;
      }
      else{
    Current[row][col] = (Previous[row-1][col] + Previous[row][col-1] + Previous[row][col] + Previous[row][col+1] + Previous[row+1][col])/(static_cast<double>(5));
      }
    }
  }
  
  for (row = 0; row < block; ++row) {
    for (col= 0; col < block; ++col) {
      Previous[row][col] = Current[row][col];
    }
  }  
}

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1;
  }

  MPI_Init(&argc, &argv);
  
  
  long N, K;
  N = atol(argv[1]);
  K = atol(argv[2]);

 int rank, size;
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Comm_size(MPI_COMM_WORLD, &size);

 long sqrt_procs = sqrt(size);
 long block = N/sqrt_procs;
 
 long row_rank = rank/sqrt_procs;
 long col_rank = rank%sqrt_procs;
  
 double** Previous = new double*[block];
 double** Current = new double*[block];

 for(long i = 0;i < block; ++i){
    Previous[i]=(double*)malloc(block * sizeof(double));
    Current[i]=(double*)malloc(block * sizeof(double));
  }


 for (long row = row_rank*block; row < (row_rank+1)*block; row++) {
   for (long col = col_rank*block; col < (col_rank+1)*block; col++) {
     Current[(row-(row_rank*block))][(col-(col_rank*block))] = generate2DHeat(N, row, col);
     Previous[(row-(row_rank*block))][(col-(col_rank*block))] = Current[(row-(row_rank*block))][(col-(col_rank*block))];
   }
 }

 MPI_Request* request;
 MPI_Status* status;
  
 double* send_left = new double[block];
 double* move_left = new double[block];
 double* send_right = new double[block];
 double* move_right = new double[block];
 double* move_up = new double[block];
 double* move_down = new double[block]; 

  MPI_Barrier(MPI_COMM_WORLD);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  for(long iter = 1; iter <= K; ++iter){

    for(long i = 0;i<block;i++){
      move_left[i]=Previous[i][0];
      move_right[i]=Previous[i][(block-1)];
      move_up[i]=Previous[0][i];
      move_down[i]=Previous[(block-1)][i];
    }
    
    if(size == 1){  
      
      calculate_2d_heat(block, Current, Previous, move_left, move_right, move_up, move_down);

      check2DHeat(Current, block, rank, sqrt_procs, iter);
      
    }
    else{
      
      if(row_rank == 0 && col_rank == 0){
    
    for(long i = 0; i < block; i++){
      send_right[i] = Previous[i][block-1];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);

    MPI_Waitall(4, request, status);
      }
      
      else if(row_rank == 0 && col_rank == (sqrt_procs-1)){
    
    for(long i = 0; i < block; i++){
      send_left[i] = Previous[i][0];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);
    
    MPI_Waitall(4, request, status);
      }
      
      else if(row_rank == (sqrt_procs-1) && col_rank == 0){
    
    for(long i = 0; i < block; i++){
      send_right[i] = Previous[i][block-1];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
    
    MPI_Waitall(4, request, status);
      }
      
      else if(row_rank == (sqrt_procs - 1) && col_rank == (sqrt_procs - 1)){
    
    for(long i = 0; i < block; i++){
      send_left[i] = Previous[i][0];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[3]);

    MPI_Waitall(4, request, status);
      }
      
      else if(row_rank == 0){
    
    for(long i = 0; i < block; i++){
      send_left[i] = Previous[i][0];
      send_right[i] = Previous[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);
    
    MPI_Waitall(6, request, status);
      }
      
      else if(row_rank == (sqrt_procs-1)){
    
    for(long i = 0; i < block; i++){
      send_left[i] = Previous[i][0];
      send_right[i] = Previous[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);

    MPI_Waitall(6, request, status);
      }
      
      
      else if(col_rank == 0){

    for(long i = 0; i < block; i++){
      send_right[i] = Previous[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];
    
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);
    MPI_Waitall(6, request, status);
      }

      else if(col_rank == (sqrt_procs-1)){
    
    for(long i = 0; i < block; i++){
      send_left[i]=Previous[i][0];
    }
    
    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[5]);
    
    MPI_Waitall(6, request, status);
      }
      
      else{
    
    for(long i = 0; i < block; i++){
      send_left[i]=Previous[i][0];
      send_right[i]=Previous[i][block-1];
    }

    MPI_Request* req_r;
    MPI_Request* req_s;
    req_r = new MPI_Request[4];
    req_s = new MPI_Request[4];
    status = new MPI_Status[4];
    
    MPI_Irecv(move_up, block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &req_r[0]);
    MPI_Irecv(move_down, block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &req_r[1]);
    MPI_Irecv(move_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_r[2]);
    MPI_Irecv(move_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_r[3]);
    MPI_Isend(Previous[0], block, MPI_DOUBLE, rank-sqrt_procs, 0, MPI_COMM_WORLD, &req_s[0]);
    MPI_Isend(Previous[block-1], block, MPI_DOUBLE, rank+sqrt_procs, 0, MPI_COMM_WORLD, &req_s[1]);
    MPI_Isend(send_left, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_s[2]);
    MPI_Isend(send_right, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_s[3]);
    MPI_Waitall(4, req_r, status);
      }
      
      calculate_2d_heat(block, Current, Previous, move_left, move_right, move_up, move_down);
      
      
      check2DHeat(Current, block, rank, sqrt_procs, iter);
      
    }
    
  }

  MPI_Finalize();
   
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  
  std::chrono::duration<double> elapsed_seconds = end-start;
  
  if(rank == 0){
    std::cerr<<elapsed_seconds.count()<<std::endl;
  }

  for(long i = 0; i< block; i++){
    delete Previous[i];
    delete Current[i];
  }
  delete[] Current;
  delete[] Previous;
  delete[] send_left;
  delete[] move_left;
  delete[] send_right;
  delete[] move_right;
  delete[] move_up;
  delete[] move_down;
  
  return 0;
}