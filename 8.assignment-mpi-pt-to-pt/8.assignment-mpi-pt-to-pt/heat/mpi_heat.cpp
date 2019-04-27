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

void calculate_2d_heat(long block, double** MainArray, double** firstArray, double* leftMArr, double* rightMArr, double* upArr, double* downArr){

  int row, col;
  
  for (row = 0; row < block; ++row) { 
    for (col = 0; col < block; ++col) {
      if(row == 0){
    MainArray[row][col] = (upArr[col] + firstArray[row][col-1] + firstArray[row][col] + firstArray[row][col+1] + firstArray[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1){
    MainArray[row][col] = (downArr[col] + firstArray[row-1][col] + firstArray[row][col-1] + firstArray[row][col] + firstArray[row][col+1])/(static_cast<double>(5));
      }
      else if(col == 0){
    MainArray[row][col] = (leftMArr[row] + firstArray[row-1][col] + firstArray[row][col] + firstArray[row][col+1] + firstArray[row+1][col])/(static_cast<double>(5));
      }
      else if(col == block-1){
    MainArray[row][col] = (rightMArr[row] + firstArray[row-1][col] + firstArray[row][col-1] + firstArray[row][col]+ firstArray[row+1][col])/(static_cast<double>(5));
      }
      else if(row == 0 && col == 0){
    MainArray[row][col] = (leftMArr[row] + upArr[col] + firstArray[row][col] + firstArray[row][col+1] + firstArray[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1 && col == 0){
    MainArray[row][col] = (leftMArr[row] + downArr[col] + firstArray[row-1][col] + firstArray[row][col] + firstArray[row][col+1])/(static_cast<double>(5));
      }
      else if(row == 0 && col == block-1){
    MainArray[row][col] = (rightMArr[row] + upArr[col] + firstArray[row][col-1] + firstArray[row][col] + firstArray[row+1][col])/(static_cast<double>(5));
      }
      else if(row == block-1 && col == block-1){
    MainArray[row][col] = (rightMArr[row] + downArr[col] + firstArray[row-1][col] + firstArray[row][col-1] + firstArray[row][col])/(static_cast<double>(5));;
      }
      else{
    MainArray[row][col] = (firstArray[row-1][col] + firstArray[row][col-1] + firstArray[row][col] + firstArray[row][col+1] + firstArray[row+1][col])/(static_cast<double>(5));
      }
    }
  }
  
  for (row = 0; row < block; ++row) {
    for (col= 0; col < block; ++col) {
      firstArray[row][col] = MainArray[row][col];
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

 long sqrtSize = sqrt(size);
 long block = N/sqrtSize;
 
 long global_ir = rank/sqrtSize;
 long global_jr = rank%sqrtSize;
  
 double** firstArray = new double*[block];
 double** MainArray = new double*[block];

 for(long i = 0;i < block; ++i){
    firstArray[i]=(double*)malloc(block * sizeof(double));
    MainArray[i]=(double*)malloc(block * sizeof(double));
  }


 for (long global_i = global_ir*block; global_i < (global_ir+1)*block; global_i++) {
   for (long global_j = global_jr*block; global_j < (global_jr+1)*block; global_j++) {
     MainArray[(global_i-(global_ir*block))][(global_j-(global_jr*block))] = generate2DHeat(N, global_i, global_j);
     firstArray[(global_i-(global_ir*block))][(global_j-(global_jr*block))] = MainArray[(global_i-(global_ir*block))][(global_j-(global_jr*block))];
   }
 }

 MPI_Request* request;
 MPI_Status* status;
  
 double* leftSArr = new double[block];
 double* leftMArr = new double[block];
 double* rightSArr = new double[block];
 double* rightMArr = new double[block];
 double* upArr = new double[block];
 double* downArr = new double[block]; 

  MPI_Barrier(MPI_COMM_WORLD);
  
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

  for(long iter = 1; iter <= K; ++iter){

    for(long i = 0;i<block;i++){
      leftMArr[i]=firstArray[i][0];
      rightMArr[i]=firstArray[i][(block-1)];
      upArr[i]=firstArray[0][i];
      downArr[i]=firstArray[(block-1)][i];
    }
    
    if(size == 1){  
      
      calculate_2d_heat(block, MainArray, firstArray, leftMArr, rightMArr, upArr, downArr);

      check2DHeat(MainArray, block, rank, sqrtSize, iter);
      
    }
    else{
      
      if(global_ir == 0 && global_jr == 0){
    
    for(long i = 0; i < block; i++){
      rightSArr[i] = firstArray[i][block-1];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[3]);

    MPI_Waitall(4, request, status);
      }
      
      else if(global_ir == 0 && global_jr == (sqrtSize-1)){
    
    for(long i = 0; i < block; i++){
      leftSArr[i] = firstArray[i][0];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[3]);
    
    MPI_Waitall(4, request, status);
      }
      
      else if(global_ir == (sqrtSize-1) && global_jr == 0){
    
    for(long i = 0; i < block; i++){
      rightSArr[i] = firstArray[i][block-1];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[3]);
    
    MPI_Waitall(4, request, status);
      }
      
      else if(global_ir == (sqrtSize - 1) && global_jr == (sqrtSize - 1)){
    
    for(long i = 0; i < block; i++){
      leftSArr[i] = firstArray[i][0];
    }

    request = new MPI_Request[4];
    status = new MPI_Status[4];

    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[3]);

    MPI_Waitall(4, request, status);
      }
      
      else if(global_ir == 0){
    
    for(long i = 0; i < block; i++){
      leftSArr[i] = firstArray[i][0];
      rightSArr[i] = firstArray[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);
    
    MPI_Waitall(6, request, status);
      }
      
      else if(global_ir == (sqrtSize-1)){
    
    for(long i = 0; i < block; i++){
      leftSArr[i] = firstArray[i][0];
      rightSArr[i] = firstArray[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);

    MPI_Waitall(6, request, status);
      }
      
      
      else if(global_jr == 0){

    for(long i = 0; i < block; i++){
      rightSArr[i] = firstArray[i][block-1]; 
    }

    request = new MPI_Request[6];
    status = new MPI_Status[6];
    
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &request[5]);
    MPI_Waitall(6, request, status);
      }

      else if(global_jr == (sqrtSize-1)){
    
    for(long i = 0; i < block; i++){
      leftSArr[i]=firstArray[i][0];
    }
    
    request = new MPI_Request[6];
    status = new MPI_Status[6];

    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[1]);
    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &request[3]);
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &request[4]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &request[5]);
    
    MPI_Waitall(6, request, status);
      }
      
      else{
    
    for(long i = 0; i < block; i++){
      leftSArr[i]=firstArray[i][0];
      rightSArr[i]=firstArray[i][block-1];
    }

    MPI_Request* req_r;
    MPI_Request* req_s;
    req_r = new MPI_Request[4];
    req_s = new MPI_Request[4];
    status = new MPI_Status[4];
    
    MPI_Irecv(upArr, block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &req_r[0]);
    MPI_Irecv(downArr, block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &req_r[1]);
    MPI_Irecv(leftMArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_r[2]);
    MPI_Irecv(rightMArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_r[3]);
    MPI_Isend(firstArray[0], block, MPI_DOUBLE, rank-sqrtSize, 0, MPI_COMM_WORLD, &req_s[0]);
    MPI_Isend(firstArray[block-1], block, MPI_DOUBLE, rank+sqrtSize, 0, MPI_COMM_WORLD, &req_s[1]);
    MPI_Isend(leftSArr, block, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_s[2]);
    MPI_Isend(rightSArr, block, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_s[3]);
    MPI_Waitall(4, req_r, status);
      }
      
      calculate_2d_heat(block, MainArray, firstArray, leftMArr, rightMArr, upArr, downArr);
      
      
      check2DHeat(MainArray, block, rank, sqrtSize, iter);
      
    }
    
  }

  MPI_Finalize();
   
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  
  std::chrono::duration<double> elapsed_seconds = end-start;
  
  if(rank == 0){
    std::cerr<<elapsed_seconds.count()<<std::endl;
  }

  for(long i = 0; i< block; i++){
    delete firstArray[i];
    delete MainArray[i];
  }
  delete[] MainArray;
  delete[] firstArray;
  delete[] leftSArr;
  delete[] leftMArr;
  delete[] rightSArr;
  delete[] rightMArr;
  delete[] upArr;
  delete[] downArr;
  
  return 0;
}