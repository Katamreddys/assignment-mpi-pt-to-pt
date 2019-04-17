#include <mpi.h>
#include <iostream>
#include <stdio.h>

int main (int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <value>"<<std::endl;
    return -1;
  }
 int a = atoi(argv[1]);
 int rank;
 MPI_Init (&argc, &argv);
 MPI_Comm_rank (MPI_COMM_WORLD, &rank);
 if(rank==0){
 	MPI_Send(&a,1,MPI_INT,1,0,MPI_COMM_WORLD);
        MPI_Recv(&a,1,MPI_INT,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	std::cout<<a<<std::endl;
 }
else
{
	MPI_Recv(&a,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	a = a+2;
	MPI_Send(&a,1,MPI_INT,0,0,MPI_COMM_WORLD);
}
  
	MPI_Finalize();

  return 0;
}