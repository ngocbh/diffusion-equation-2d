/*
 *  Problem : main.c
 *  Description : 
 *  Created by ngocjr7 on 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>


using namespace std;
const long N = 100000 + 7;
const long INF = 1000000000 + 7;
const long MODULE = 1000000000 + 7;

#define m 1000
#define n 1000
#define epsilon 0.001
#define tolerance 0.001

void copy_array(float* C,float* Cn,int M,int N)
{
	for (int i = 0; i < M+1; i++)
		for (int j = 0; j < N; j++)
			*(C+i*n+j) = *(Cn+i*n+j);
}

void print_array(float* C)
{
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%.1f\t", *(C+i*n+j));
		printf("\n");
	}
}

void KhoiTao(float *C)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if ( i >= (m/2-5) && i < (m/2+5) && j >= (n/2-5) && j < (n/2+5) ) 
				*(C+i*n+j) = 80.0;
			else 
				*(C+i*n+j) = 25.0; 
}

void calc_red_back(float *C,float *Cn,int M,int N,int rb,int rank,float* delta)
{
	for (int i = 1; i <= M; i++)
		for (int j = 0; j < N; j++) {
			int ii = M*rank + i - 1;
			int ij = j;
			if ( (ii+ij)%2 == rb ) {
				if (ii >= (m/2-5) && ii < (m/2+5) && ij >= (n/2-5) && ij < (n/2+5)) {
					*(Cn+i*N+j) = 80.0;
				}
				else if ( ii == 0 || ii == m-1 || ij == 0 || ij == n-1 ) {
					*(Cn+i*N+j) = 25.0;
				}
				else {
					*(Cn+i*N+j) = 0.25 * (*(C+(i-1)*N+j) + *(C+(i+1)*N+j) + *(C+i*N+j-1) + *(C+i*N+j+1));
				}

				if ( abs(*(Cn+i*N+j) - *(C+i*N+j)) > tolerance ) 
					*delta = abs(*(Cn+i*N+j) - *(C+i*N+j));
			} else {
				*(Cn+i*N+j) = *(C+i*N+j);
			}
		}
}

int main(int argc,char *argv[])
{
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout);

	int NP, rank, mc;
	float *C;
	C = (float*)malloc(m*n*sizeof(float));
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&NP);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank == 0) {
    	KhoiTao(C);
    	// print_array(C);
    }

    float delta;
	float *Cn, *Cs, *Cl, *Cr;
	MPI_Status status;

	mc = m/NP;
	Cn = (float*)malloc((mc+2)*n*sizeof(float));
	Cs = (float*)malloc((mc+2)*n*sizeof(float));

	MPI_Scatter(C,mc*n,MPI_FLOAT,Cs+n,mc*n,MPI_FLOAT,0,MPI_COMM_WORLD);

	do {
		delta = 0.0;
		//Calc red node
		if (rank == 0) {
            MPI_Send(Cs+mc*n,n,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
        } else if ( rank == NP - 1) {
            MPI_Recv(Cs,n,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
        } else {
            MPI_Send(Cs+mc*n,n,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
            MPI_Recv(Cs,n,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
        }
        if (rank == NP-1) {
            MPI_Send(Cs+n,n,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
        } else if ( rank == 0 ) {
            MPI_Recv(Cs+(mc+1)*n,n,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
        } else {
            MPI_Send(Cs+n,n,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
            MPI_Recv(Cs+(mc+1)*n,n,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
        }
		calc_red_back(Cs,Cn,mc,n,0,rank,&delta);
		copy_array(Cs,Cn,mc,n);
		//Calc back node
		if (rank == 0) {
            MPI_Send(Cs+mc*n,n,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
        } else if ( rank == NP - 1) {
            MPI_Recv(Cs,n,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
        } else {
            MPI_Send(Cs+mc*n,n,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
            MPI_Recv(Cs,n,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
        }
        if (rank == NP-1) {
            MPI_Send(Cs+n,n,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
        } else if ( rank == 0 ) {
            MPI_Recv(Cs+(mc+1)*n,n,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
        } else {
            MPI_Send(Cs+n,n,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
            MPI_Recv(Cs+(mc+1)*n,n,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
        }
		calc_red_back(Cs,Cn,mc,n,1,rank,&delta);
		copy_array(Cs,Cn,mc,n);

		//get delta
		float new_delta = 0;
		MPI_Reduce(&delta,&new_delta,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
		if ( rank == 0 ) {
			delta = new_delta;
		}
		MPI_Bcast(&delta,1,MPI_FLOAT,0,MPI_COMM_WORLD);
	} while (delta > tolerance); 


	MPI_Gather(Cs+n,mc*n,MPI_FLOAT,C,mc*n,MPI_FLOAT,0,MPI_COMM_WORLD);

	if (rank == 0) {
		print_array(C);
	}

	MPI_Finalize();

	return 0;
}
