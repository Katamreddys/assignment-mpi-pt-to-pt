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


int rank, numprocs;

MPI_Init(&argc, &argv);

MPI_Status status;

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

MPI_Comm_size(MPI_COMM_WORLD, &numprocs);




double x,integration;

int fid = atoi(argv[1]);

int a = atoi(argv[2]);

int b = atoi(argv[3]);

int n = atoi(argv[4]);

int intensity = atoi(argv[5]);

int first = 0, last = 0, end_proc = -1;

auto start = chrono::high_resolution_clock::now();


int count = numprocs;

int granularity = (n/numprocs)/8;


if(numprocs > 1){

//for the Master code

if(rank == 0){

for(int i = 1; i < numprocs; ++i)
{

int j = i;

	if(n >= first + granularity)
	{

		MPI_Send(&first, 1, MPI_INT, i, j, MPI_COMM_WORLD);

		first += granularity;

	}

	else{

		if(count == numprocs)
		{

			MPI_Send(&first, 1, MPI_INT, i, j, MPI_COMM_WORLD);

			count--;

		}

		else{

			count--;

			MPI_Send(&end_proc, 1, MPI_INT, i, j, MPI_COMM_WORLD);

		    }

            }

}

while(count != 1)
{

double sum = 0.0;

MPI_Recv(&sum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

integration += sum;

	if(first < n)
	{

		MPI_Send(&first, 1, MPI_INT, status.MPI_SOURCE, rank, MPI_COMM_WORLD);

		first += granularity;

	}

	else
	{

		count--;

		MPI_Send(&end_proc, 1, MPI_INT, status.MPI_SOURCE, rank, MPI_COMM_WORLD);

	}

}

integration = integration * ((float)(b-a)/n);

}

// for the Worker code
	else
	{

		while(1)
		{

			double sum = 0.0;

			MPI_Recv(&first, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if(first == -1)

				break;

			last = first + granularity;

			if(last > n)

				last = n;

	for(int i = first; i < last; ++i)
	{

		x = (a + (i + 0.5) * ((float)(b-a)/n));

		switch(fid){

			case 1:

				sum += f1(x, intensity);
	
				break;

			case 2:

				sum += f2(x, intensity);

				break;

			case 3:

				sum += f3(x, intensity);

				break;

			case 4:

				sum += f4(x, intensity);

				break;

		        default: exit;

			    }

	}

MPI_Send(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

}

}

}


else

{

double sum = 0.0;

	for(int i = 0; i < n; ++i){

		x = (a + (i + 0.5) * ((float)(b-a)/n));

		switch(fid){

			case 1:

				sum += f1(x, intensity);

				break;

			case 2:

				sum += f2(x, intensity);

				break;

			case 3:

				sum += f3(x, intensity);
	
				break;

			case 4:

				sum += f4(x, intensity);

				break;

			default: exit;

				}

	}

integration = sum * ((float)(b-a)/n);

}

auto end = chrono::high_resolution_clock::now() - start;

if(rank == 0){

cout << integration << endl;

cerr << chrono::duration<double>(end).count()<< endl;

}

MPI_Finalize();

return 0;

}