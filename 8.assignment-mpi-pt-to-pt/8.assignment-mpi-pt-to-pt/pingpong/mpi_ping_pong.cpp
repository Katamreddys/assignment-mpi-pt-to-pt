#include <mpi.h>
#include <iostream>
using namespace std;

int main (int argc, char* argv[]) {

  if (argc < 2) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <element>"<<std::endl;
    return -1;
  }
  int element;
  MPI_Init(&argc,&argv);
  element = stoi(argv[1]);

  // rank and size
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int number;
  if (rank == 0) {
      //element
    number = element;
      //send and receive
    MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    cout<<"Sending from 0 to 1"<<endl;
    MPI_Recv(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"Received after adding"<<endl;
    cout<<"After adding 2: " <<number<<endl;
  } else if (rank == 1) {
      
      //receive add and send
    MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    cout<<"Received to 1 from 0"<<endl;

    number = number + 2;
    MPI_Send(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    cout<<"Sending from 1 to 0 after adding 2"<<endl;

  }


  MPI_Finalize();
  return 0;
}
