/*
 *  Problem : main.c
 *  Description : 
 *  Created by ngocjr7 on 
*/

#include <stdlib.h>
#include <stdio.h>

using namespace std;
const long N = 100000 + 7;
const long INF = 1000000000 + 7;
const long MODULE = 1000000000 + 7;
typedef pair<int,int> ii;

#define m 20
#define n 20
#define epsilon 0.001
#define tolerance 0.001

void KhoiTao(float *C)
{
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if ( i >= (m/2-5) && i < (m/2+5) && j >= (n/2-5) && j < (n/2+5) ) 
				*(C+i*n+j) = 80.0;
			else 
				*(C+i*n+j) = 25.0; 
}

int main()
{
	freopen("in.txt","r",stdin);

	

	return 0;
}
