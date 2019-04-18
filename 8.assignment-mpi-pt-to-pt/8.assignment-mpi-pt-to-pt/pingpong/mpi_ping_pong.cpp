#include <mpi.h>

// C++ std lib
#include <iostream>

const int MASTER = 0;

int main(int argc, char *argv[])
{
    MPI :: Init(argc, argv);
    MPI :: COMM_WORLD.Set_errhandler(MPI :: ERRORS_THROW_EXCEPTIONS);

    try {

        int mes = 0;
        int rank = MPI :: COMM_WORLD.Get_rank();
        int n_rank = rank + 1;
        int p_rank = rank - 1;
        int m_rank = MPI :: COMM_WORLD.Get_size() - 1;

        if (rank == MASTER) {
            std :: cerr << "I am " << rank << std :: endl;

            if (rank != m_rank)
                while (n_rank <= m_rank) {

                    std :: cerr << "send from " << rank << "| to " << n_rank << std :: endl;
                    MPI :: COMM_WORLD.Send(&n_rank, 1, MPI :: INT, n_rank, 1);
                    MPI :: COMM_WORLD.Recv(&mes, 1, MPI :: INT, n_rank, 1);
                    std :: cerr << "recv from " << n_rank << "| to " << rank << std :: endl;
                    ++n_rank;
                }
        } else {

            MPI :: COMM_WORLD.Recv(&mes, 1, MPI :: INT, MASTER, 1);
            std :: cerr << "I am " << rank << "| my rank from master =  " << mes <<  std :: endl;
            std :: cerr << "send to MASTER"  << "| I'm " << rank << std :: endl;
            MPI :: COMM_WORLD.Send(&rank, 1, MPI :: INT, MASTER, 1);
        }
    } catch (MPI :: Exception e) {

        std::cout << "MPI ERROR: " << e.Get_error_code() \
        << " - " << e.Get_error_string() << std :: endl;
    }

    MPI :: Finalize();
    return 0;
}