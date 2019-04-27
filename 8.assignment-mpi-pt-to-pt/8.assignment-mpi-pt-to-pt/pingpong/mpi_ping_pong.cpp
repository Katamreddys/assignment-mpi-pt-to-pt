#include <mpi.h>
#include <iostream>
using namespace std;

int main (int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <element>"<<std::endl;
    return -1;
  }
  int element,value,mpi_rank, mpi_size;
  MPI_Init(&argc,&argv);
  element = stoi(argv[1]);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank == 0) {
    value = element;
    MPI_Send(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    cout<<"0 Send to 1"<<endl;
    MPI_Recv(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"Recive After addind"<<endl;
    cout<<"Final Value: " <<value<<endl;
  } else if (mpi_rank == 1) {
    MPI_Recv(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"Received by 1 from 0"<<endl;
    value = value + 2;//add 2
    MPI_Send(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    cout<<"Sending from 1 to 0 after adding 2"<<endl;
  }
  MPI_Finalize();
  return 0;
}
