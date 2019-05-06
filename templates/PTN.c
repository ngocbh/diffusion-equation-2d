/*
	Created by cosenkid on 08042019
*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define M       20
#define Time    1
#define dt      0.01
#define dx      0.1
#define D       0.1

void KhoiTao(float *T)
{
        int i;
        for (i = 0; i < M; i++)
        {
                *(T+i) = 25.0;
        }
}

void DHBH(float *Ts,float *dTs,float *Tl,float *Tr,int Mc)
{
        int i;
        float c, l, r;
        for (i = 0; i < Mc; i++)
        {
                c = *(Ts+i);
                l = (i == 0)? *Tl : *(Ts + i - 1);
                r = (i == Mc-1)? *Tr : *(Ts + i + 1);
                *(dTs + i) = D*(l-2*c+r) / (dx * dx);
        }
}

int main(int argc,char *argv[])
{
        int i;
        int NP, rank, Mc;
        MPI_Status status;
        float *T, *dT, *Ts, *dTs;
        float *Tl, *Tr;
        T = (float*)malloc(M*sizeof(float));
        dT = (float*)malloc(M*sizeof(float));
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&NP);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0) {
                KhoiTao(T);     
        }
        Mc = M/NP;
        Ts = (float*)malloc(Mc*sizeof(float));
        dTs = (float*)malloc(Mc*sizeof(float));

	Tl = (float*)malloc(sizeof(float));
	Tr = (float*)malloc(sizeof(float));

        for (int i = 0; i < M; i++) 
                *(T+i) = i;

        MPI_Scatter(T,Mc,MPI_FLOAT,Ts,Mc,MPI_FLOAT,0,MPI_COMM_WORLD);
        float t = 0 ;

        printf("rank = %d\n", rank);
        for (int i = 0; i < Mc; i++) {
                printf("[%d %f] ",rank,*(Ts+i));
        }
        fflush(stdout);


        while ( t < Time ) {
                // printf("%d\n",rank);
                if (rank == 0) {
                        *Tl = 100.0;
                        MPI_Send(Ts+Mc-1,1,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
                } else if ( rank == NP - 1) {
                        MPI_Recv(Tl,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
                } else {
                        MPI_Send(Ts+Mc-1,1,MPI_FLOAT,rank+1,rank,MPI_COMM_WORLD);
                        MPI_Recv(Tl,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD,&status);
                }
                if (rank == NP-1) {
                        *Tr = 25.0;
                        MPI_Send(Ts,1,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
                } else if ( rank == 0 ) {
                        MPI_Recv(Tr,1,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
                } else {
                        MPI_Send(Ts,1,MPI_FLOAT,rank-1,rank,MPI_COMM_WORLD);
                        MPI_Recv(Tr,1,MPI_FLOAT,rank+1,rank+1,MPI_COMM_WORLD,&status);
                }
                DHBH(Ts,dTs,Tl,Tr,Mc);

                for (i =  0; i < Mc; i++) {
                        *(Ts+i) = *(Ts+i) + *(dTs+i) * dt;
                }
                // printf("hala one: %d\n",rank);
                t = t + dt;
        }

        MPI_Gather(Ts,Mc,MPI_FLOAT,T,Mc,MPI_FLOAT,0,MPI_COMM_WORLD);
        if (rank == 0) {
                for (i = 0; i < M; i++)
                        printf("%5.2f ",*(T+i));
        }
        MPI_Finalize();
        return 0;
}



