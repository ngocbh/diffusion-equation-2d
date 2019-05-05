#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#define MAX 1000000

int main(int argc, char** argv) {
	int rank, size;
	MPI_Status status;
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);
	int s[MAX];
	int i;
	for (i = 0; i < MAX; i++) 
		s[i] = i;
	int r[MAX];
	if (rank == 0) {
	MPI_Send(s,MAX,MPI_INT,1-rank, rank, MPI_COMM_WORLD);
	MPI_Recv(r,MAX,MPI_INT,1-rank, 1 - rank, MPI_COMM_WORLD, &status);
	} else {
	        MPI_Recv(r,MAX,MPI_INT,1-rank, 1 - rank, MPI_COMM_WORLD, &status);
		MPI_Send(s,MAX,MPI_INT,1-rank, rank, MPI_COMM_WORLD);
	}
	printf ("Process %d received: ", rank);
	for (i = 0; i < MAX; i++) {
		printf ("%d ", r[i]);
	}
	printf("\n");
	MPI_Finalize();
	return 0;
}
