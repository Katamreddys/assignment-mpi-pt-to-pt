#include <mpi.h>

// C++ std lib
#include <iostream>

const int MASTER = 0;

int main(int argc, char *argv[])
{
    MPI :: Init(argc, argv);
    MPI :: COMM_WORLD.Set_errhandler(MPI :: ERRORS_THROW_EXCEPTIONS);

    try {

        int message = 0;
        int mpi_rank = MPI :: COMM_WORLD.Get_rank();
        int rankN = mpi_rank + 1;
        int rankM = MPI :: COMM_WORLD.Get_size() - 1;

        if (mpi_rank == MASTER) {
            std :: cerr << "I am " << mpi_rank << std :: endl;

            if (mpi_rank != rankM)
                while (rankN <= rankM) {

                    std :: cerr << "send from " << mpi_rank << "| to " << rankN << std :: endl;
                    MPI :: COMM_WORLD.Send(&rankN, 1, MPI :: INT, rankN, 1);
                    MPI :: COMM_WORLD.Recv(&message, 1, MPI :: INT, rankN, 1);
                    std :: cerr << "recv from " << rankN << "| to " << mpi_rank << std :: endl;
                    ++rankN;
                }
        } else {

            MPI :: COMM_WORLD.Recv(&message, 1, MPI :: INT, MASTER, 1);
            std :: cerr << "I am " << mpi_rank << "| my mpi_rank from master =  " << message <<  std :: endl;
            std :: cerr << "send to MASTER"  << "| I'm " << mpi_rank << std :: endl;
            MPI :: COMM_WORLD.Send(&mpi_rank, 1, MPI :: INT, MASTER, 1);
        }
    } catch (MPI :: Exception e) {

        std::cout << "MPI ERROR: " << e.Get_error_code() \
        << " - " << e.Get_error_string() << std :: endl;
    }

    MPI :: Finalize();
    return 0;
}