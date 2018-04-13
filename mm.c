/* mm.c 
 * Sample matrix multiplication in C 
 * Author: Philip Papadopoulos
 * Email: ppapadopoulos@ucsd.edu
 * UCSD Course: cse160 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
double cpu_time ( void );
void mult(double **, double **, double **, int, int, int);
void printMatrices(double **, double **, double **, int, int, int);
void usage( void );

int main(int argc, char * argv[]) {
	if (argc != 5) usage();

	int N = atoi(argv[1]);
	int M = atoi(argv[2]);
	int P = atoi(argv[3]); 
	int rounds = atoi(argv[4]);

	double **A, **B, **C;
	double *Alinear, *Blinear, *Clinear;
	A = (double **) calloc(N, sizeof (double **));
	B = (double **) calloc(M, sizeof (double **));
	C = (double **) calloc(N, sizeof (double **));
	Alinear = calloc(N * M, sizeof (double));
	Blinear = calloc(M * P, sizeof (double));
	Clinear = calloc(N * P, sizeof (double));

	int i;
	/* populate A, B ,C so that we can use  A[i][j] addressing */
	for (i = 0; i < N ; i++) { 
		A[i] = Alinear + i * M;
		C[i] = Clinear + i * P;
	}

	for (i = 0; i < M; i++) {
		B[i] = Blinear + i * P;
	}

	for (i = 0 ; i < N * M; i++) {
		Alinear[i] = drand48();
	}

	for (i = 0; i < M * P; i++) {
		Blinear[i] = drand48();
	}


	double startTime, endTime;
	printf("Multiplying (%d x %d) and (%d x %d) random matrices\n", N, M, M, P);
	startTime = cpu_time();
    
	int row, col, k;
	double sum;
	for (i = 0; i < rounds; i ++) {
		mult(A, B, C, N, M, P);
	}
    
	endTime = cpu_time(); 
	double avgTime = (endTime - startTime) / (double) rounds;
	printf("Multiplied (%d x %d) and (%d x %d) in %f seconds (%d trials)\n\n", N, M, M, P, avgTime, rounds);

	// printMatrices(A, B, C, N, M, P);

	free(A);
	free(Alinear);
	free(B);
	free(Blinear);
	free(C);
	free(Clinear);
	return 0;
}

/** printMatrices
 * Print contents of matrices A, B, and C to stdout.
 * Perform after multiplying matrices, C = A x B.
*/
void printMatrices(double **A, double **B, double**C, int N, int M, int P) {
	int i, j;
	printf("Matrix A contents:\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			printf("%6.4f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
	printf("Matrix B contents:\n");
	for (i = 0; i < M; i++) {
                for (j = 0; j < P; j++) {
                        printf("%6.4f ", B[i][j]);
                }
                printf("\n");
        }
	printf("\n\n");
	printf("Matrix C contents:\n");
	for (i = 0; i < N; i++) {
                for (j = 0; j < P; j++) {
                        printf("%6.4f ", C[i][j]);
                }
                printf("\n");
        }
	printf("\n\n"); 
}

/** mult
 * Multiply matrices, C = A x B.
 * Assumes Matrices are already allocated.
 * values are not checked for sanity
*/
void mult(double **A, double **B, double **C, int N, int M, int P) {
	int k, row, col;
	double sum;
	for (row = 0; row < N; row++) {
		for (col = 0; col < P; col++) {
			sum = 0.0;
			for (k = 0 ; k < M; k++) {
				sum += A[row][k] * B[k][col];
			}

			C[row][col] = sum;
		}
	}
}

/** usage
 * Print usage message to stderr
*/
void usage(void) {
	fprintf (stderr, "usage: mm <N> <M> <P> <trials>\n");
	exit(-1);
}

/** cpu_time
 * Purpose:
 * CPU_TIME returns the current reading on the CPU clock.
 * Licensing:
 * This code is distributed under the GNU LGPL license.
 * Modified:
 * 06 June 2005
 * Author:
 * John Burkardt
 * Parameters:
 * Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
double cpu_time ( void ) {
	double value;
	value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;
	return value;
}
/* vim: set ts=4
*/
