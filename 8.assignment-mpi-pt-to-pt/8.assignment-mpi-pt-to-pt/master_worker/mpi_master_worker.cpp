#include <mpi.h>

#include <iostream>

#include <chrono>


using namespace std;


#ifdef __cplusplus

extern "C" {

#endif


float f1(float x, int intensity);

float f2(float x, int intensity);

float f3(float x, int intensity);

float f4(float x, int intensity);


#ifdef __cplusplus

}

#endif


//int granularity = 10000;


int main (int argc, char* argv[]) {


if (argc < 6) {

std::cerr<<"usage: mpirun "<<argv[0]<<" <functionid> <a> <b> <n> <intensity>"<<std::endl;

return -1;

}


int mpi_rank, proc_count;

MPI_Init(&argc, &argv);

MPI_Status mpi_status;

MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

MPI_Comm_size(MPI_COMM_WORLD, &proc_count);




double x,total_sum;

int function_id = atoi(argv[1]);

int a = atoi(argv[2]);

int b = atoi(argv[3]);

int n = atoi(argv[4]);

int intensity = atoi(argv[5]);

int begin = 0, stop = 0, proc_stop = -1;

auto time_start = chrono::high_resolution_clock::now();


int message_count = proc_count;

int granularity = (n/proc_count)/8;


if(proc_count > 1){

//for the Master code

if(mpi_rank == 0){

for(int i = 1; i < proc_count; ++i)
{

int j = i;

	if(n >= begin + granularity)
	{

		MPI_Send(&begin, 1, MPI_INT, i, j, MPI_COMM_WORLD);

		begin += granularity;

	}

	else{

		if(message_count == proc_count)
		{

			MPI_Send(&begin, 1, MPI_INT, i, j, MPI_COMM_WORLD);

			message_count--;

		}

		else{

			message_count--;

			MPI_Send(&proc_stop, 1, MPI_INT, i, j, MPI_COMM_WORLD);

		    }

            }

}

while(message_count != 1)
{

double partial_sum = 0.0;

MPI_Recv(&partial_sum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);

total_sum += partial_sum;

	if(begin < n)
	{

		MPI_Send(&begin, 1, MPI_INT, mpi_status.MPI_SOURCE, mpi_rank, MPI_COMM_WORLD);

		begin += granularity;

	}

	else
	{

		message_count--;

		MPI_Send(&proc_stop, 1, MPI_INT, mpi_status.MPI_SOURCE, mpi_rank, MPI_COMM_WORLD);

	}

}

total_sum = total_sum * ((float)(b-a)/n);

}

// for the Worker code
	else
	{

		while(1)
		{

			double partial_sum = 0.0;

			MPI_Recv(&begin, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);

			if(begin == -1)

				break;

			stop = begin + granularity;

			if(stop > n)

				stop = n;

	for(int i = begin; i < stop; ++i)
	{

		x = (a + (i + 0.5) * ((float)(b-a)/n));

		switch(function_id){

			case 1:

				partial_sum += f1(x, intensity);
	
				break;

			case 2:

				partial_sum += f2(x, intensity);

				break;

			case 3:

				partial_sum += f3(x, intensity);

				break;

			case 4:

				partial_sum += f4(x, intensity);

				break;

		        default: exit;

			    }

	}

MPI_Send(&partial_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

}

}

}


else

{

double partial_sum = 0.0;

	for(int i = 0; i < n; ++i){

		x = (a + (i + 0.5) * ((float)(b-a)/n));

		switch(function_id){

			case 1:

				partial_sum += f1(x, intensity);

				break;

			case 2:

				partial_sum += f2(x, intensity);

				break;

			case 3:

				partial_sum += f3(x, intensity);
	
				break;

			case 4:

				partial_sum += f4(x, intensity);

				break;

			default: exit;

				}

	}

total_sum = partial_sum * ((float)(b-a)/n);

}

auto time_end = chrono::high_resolution_clock::now() - time_start;

if(mpi_rank == 0){

cout << total_sum << endl;

cerr << chrono::duration<double>(time_end).count()<< endl;

}

MPI_Finalize();

return 0;

}