/*
 *  Problem : main.c
 *  Description :
 *  Created by ngocjr7 on
*/

#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
const long N = 100000 + 7;
const long INF = 1000000000 + 7;
const long MODULE = 1000000000 + 7;

#define m 100
#define n 100
#define epsilon 0.001
#define tolerance 0.001

void copy_array(float* C,float* Cn)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
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

void Gauss_Seidel(float *C)
{

	float delta = 0;
	float *Cn;
	Cn = (float*)malloc(n*m*sizeof(float));
    int cnt = 0, gap = 0, loop = 0, printed = 0;

	do {
		delta = 0;
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if ( (i+j)%2 == 0 ) {
					if (i >= (m/2-5) && i < (m/2+5) && j >= (n/2-5) && j < (n/2+5)) {
						*(Cn+i*n+j) = 80.0;
					}
					else if ( i == 0 || i == m-1 || j == 0 || j == n-1 ) {
						*(Cn+i*n+j) = 25.0;
					}
					else {
						*(Cn+i*n+j) = 0.25 * (*(C+(i-1)*n+j) + *(C+(i+1)*n+j) + *(C+i*n+j-1) + *(C+i*n+j+1));
					}

					if ( abs(*(Cn+i*n+j) - *(C+i*n+j)) > tolerance )
						delta = abs(*(Cn+i*n+j) - *(C+i*n+j));
				} else {
					*(Cn+i*n+j) = *(C+i*n+j);
				}

		copy_array(C,Cn);

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				if ( (i+j)%2 == 1 ) {
					if (i >= (m/2-5) && i < (m/2+5) && j >= (n/2-5) && j < (n/2+5)) {
						*(Cn+i*n+j) = 80.0;
					}
					else if ( i == 0 || i == m-1 || j == 0 || j == n-1 ) {
						*(Cn+i*n+j) = 25.0;
					}
					else {
						*(Cn+i*n+j) = 0.25 * (*(C+(i-1)*n+j) + *(C+(i+1)*n+j) + *(C+i*n+j-1) + *(C+i*n+j+1));
					}

					if ( abs(*(Cn+i*n+j) - *(C+i*n+j)) > tolerance )
						delta = abs(*(Cn+i*n+j) - *(C+i*n+j));
				} else {
					*(Cn+i*n+j) = *(C+i*n+j);
				}

		copy_array(C,Cn);
        if ( printed < 50 ) {
            print_array(C);
        } else {
			if (printed % 20 == 0) {
				print_array(C);
			}
		}
		++printed;
		// cout << "-----------------------------" << endl << delta << endl << "---------------------------" << endl;
	} while (delta > tolerance);

	print_array(C);
	cout << printed;
}

int main()
{
	freopen("in.txt","r",stdin);
	freopen("out.txt","w",stdout);
	float *C;
	C = (float*)malloc(m*n*sizeof(float));
	KhoiTao(C);
	print_array(C);
	Gauss_Seidel(C);

	return 0;
}
