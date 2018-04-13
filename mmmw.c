#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

void usage(void);
double cpu_time(void);
void printMatrices(double **, double**, double **, int, int, int);
double *generateLinearFromTwoBlocks(int, int, int, int, double **, double**);
double *generateLinearFromBlock(int, double **);
double **generateBlockFromLinear(int, double *, double **);
double **mult(int, int, double *, double **, double **, double **);

int main (int argc, char * argv[]) {
	if (argc != 6) usage();

	int busy = 1;
	int cont;
	int free = -1;
	int my_id;
	int nproc;
	int count = 0;
	int blocksize = atoi(argv[1]);
	int N = atoi(argv[2]);
	int M = atoi(argv[3]);
	int P = atoi(argv[4]);
	int currBlock = 0;
	char *output_file = argv[5];
	double sum;
	double ctime;
	double ctime1;
	double ctime2;
	double **A = (double **) calloc(N * blocksize, sizeof(double **));
	double **B = (double **) calloc(M * blocksize, sizeof(double **));
	double **C = (double **) calloc(N * blocksize, sizeof(double **));
	double **blockA = (double **) calloc(blocksize, sizeof(double **));
	double **blockB = (double **) calloc(blocksize, sizeof(double **));
	double **blockC = (double **) calloc(blocksize, sizeof(double **));
	double **blockSum = (double **) calloc(blocksize, sizeof(double **));
	double *Alinear = calloc(N * M * blocksize * blocksize, sizeof(double));
	double *Blinear = calloc(M * P * blocksize * blocksize, sizeof(double));
	double *Clinear = calloc(N * P * blocksize * blocksize, sizeof(double));
	double *blockALinear = calloc(blocksize * blocksize, sizeof(double));
	double *blockBLinear = calloc(blocksize * blocksize, sizeof(double));
	double *blockCLinear = calloc(blocksize * blocksize, sizeof(double));
	double *blockSumLinear = calloc(blocksize * blocksize, sizeof(double));
	double *blockTemp = calloc(blocksize * blocksize, sizeof(double));
	double *blockLinear = calloc(blocksize * blocksize * 2, sizeof(double));
	FILE *fp;

	MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Status status[nproc];
	MPI_Status s;
	MPI_Request request[nproc];
	MPI_Request r;

	int i;
	for (i = 0; i < blocksize; i++) {
                blockA[i] = blockALinear + i * blocksize;
                blockB[i] = blockBLinear + i * blocksize;
        	blockC[i] = blockCLinear + i * blocksize;
		blockSum[i] = blockSumLinear + i * blocksize;
        }

	// Initialize arrays in process 0
	if (my_id == 0) {
		for (i = 0; i < N * blocksize; i++) { 
                	A[i] = Alinear + i * M * blocksize;
                	C[i] = Clinear + i * P * blocksize;
        	}

        	for (i = 0; i < M * blocksize; i++) {
                	B[i] = Blinear + i * P * blocksize;
        	}

        	for (i = 0 ; i < N * M * blocksize * blocksize; i++) {
                	Alinear[i] = drand48();
        	}

        	for (i = 0; i < M * P * blocksize * blocksize; i++) {
                	Blinear[i] = drand48();
        	}
	}

	ctime1 = cpu_time();

	// Start of masterworks
	if (my_id == 0) {
		// get next block
		cont = 1;
		int j, k, blockI, blockJ, nextProc;
		for (i = 0; i < N; i++) {
			for (j = 0; j < P; j++) {
				int index = 0;
				sum = 0.0;
				for (k = 0; k < M; k++) {
					MPI_Recv(&free, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &s);
					nextProc = free;

					blockLinear = generateLinearFromTwoBlocks(blocksize, i, j, k, A, B);
					int u;
					MPI_Isend(&cont, 1, MPI_INT, nextProc, 1, MPI_COMM_WORLD, &r);
					MPI_Isend(blockLinear, blocksize * blocksize * 2, MPI_DOUBLE, nextProc, 0,  MPI_COMM_WORLD, &r);
				}
				
				// gather M blocks here and sum
				int a, g, h;
				for (g = 0; g < blocksize; g++) {
                                        for (h = 0; h < blocksize; h++) {
                                        	blockSum[g][h] = 0.0;
                                        }
                                }
				
				for (a = 0; a < M; a++) {
					MPI_Recv(blockTemp, blocksize * blocksize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &s);
					blockA = generateBlockFromLinear(blocksize, blockTemp, blockA);

					int x, y;
					for (x = 0; x < blocksize; x++) {
						for (y = 0; y < blocksize; y++) {
							blockSum[x][y] += blockA[x][y];
						}
					}
				}		

				// store summed block in C
				int startI, startJ, x, y;
				x = 0;
				y = 0;
				for (startI = i * blocksize; startI < (i * blocksize) + blocksize; startI++) {
					y = 0;
					for (startJ = j * blocksize; startJ < (j * blocksize) + blocksize; startJ++) {
						C[startI][startJ] = blockSum[x][y++];
					}
					x++;
				}
			}
		}
			
		// send no more blocks message here
		cont = 0;
		for (i = 1; i < nproc; i++) {
			MPI_Send(&cont, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
		}
		
		ctime2 = cpu_time();
		ctime = ctime2 - ctime1;

	} else {
		// check if more blocks
		while (busy == 1) {
			free = my_id;
			MPI_Send(&free, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
			MPI_Recv(&cont, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &s);
			if (cont == 1) {
				// recv block struct
				MPI_Recv(blockLinear, blocksize * blocksize * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &s);
				// perform mm
				blockC = mult(my_id, blocksize, blockLinear, blockA, blockB, blockC);
				blockTemp = generateLinearFromBlock(blocksize, blockC);
				count++;				

				// send back to process 0
				MPI_Isend(blockTemp, blocksize * blocksize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &r);
				// send no busy flag to process 0
				MPI_Isend(&free, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &r);
			} else {
				busy = 0;
				ctime2 = cpu_time();
				ctime = ctime2 - ctime1;
				MPI_Send(&ctime, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD);
				MPI_Send(&count, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
			}
		}
	}

	
	// print matrices to a fiile
	if (my_id == 0) {
		double finalTime;
		int count;
		int i;
		printf("Parameters: %d, %d, %d, %d, %d\n", blocksize, N, M, P, nproc - 1);
		for (i = 1; i < nproc; i++) {
			MPI_Recv(&finalTime, 1, MPI_DOUBLE, i, 11, MPI_COMM_WORLD, &s);
			MPI_Recv(&count, 1, MPI_INT, i, 12, MPI_COMM_WORLD, &s);
			printf("Worker %d: %d matrix multiplies performed, %15.8f time working in seconds\n", i, count, finalTime);
		}
		
		printf("Total time: %15.8f\n", ctime);

		fp = fopen(output_file, "w");
		fprintf(fp, "%d,%d,%d\n", N * blocksize, M * blocksize, P * blocksize);
		
		// Print matrix A
		int j;
		for (i = 0; i < N * blocksize; i++) {
			for (j = 0; j < M * blocksize; j++) {
				if (j != ((M * blocksize) - 1)) {
					fprintf(fp, "%15.8f, ", A[i][j]);
				}
			}
			fprintf(fp, "%15.8f\n", A[i][(M * blocksize) - 1]);
		}

		// Print matrix B
		for (i = 0; i < M * blocksize; i++) {
                        for (j = 0; j < P * blocksize; j++) {
                                if (j != ((P * blocksize) - 1)) {
                                        fprintf(fp, "%15.8f, ", B[i][j]);
                                }
                        }
                        fprintf(fp, "%15.8f\n", B[i][(P * blocksize) - 1]);
                }

		// Print matrix C
		for (i = 0; i < N * blocksize; i++) {
                        for (j = 0; j < P * blocksize; j++) {
                                if (j != ((P * blocksize) - 1)) {
                                        fprintf(fp, "%15.8f, ", C[i][j]);
                                }
                        }
                        fprintf(fp, "%15.8f\n", C[i][(P * blocksize) - 1]);
		}

		fclose(fp);
		/*free(A);
        	free(Alinear);
        	free(B);	
        	free(Blinear);
        	free(C);
        	free(Clinear);*/
	}

	MPI_Finalize();		
	return 0;
}

double *generateLinearFromBlock(int blocksize, double **block) {
	double *blockLinear = calloc(blocksize * blocksize * 2, sizeof(double));
	int i, j;
	int index = 0;
	for (i = 0; i < blocksize; i++) {
		for (j = 0; j < blocksize; j++) {
			blockLinear[index++] = block[i][j];	
		}
	}

	return blockLinear;
}

double **generateBlockFromLinear(int blocksize, double *blockLinear, double **block) {
	int i, j;
	int index = 0;
	for (i = 0; i < blocksize; i++) {
		for (j = 0; j < blocksize; j++) {
			block[i][j] = blockLinear[index++];
		}
	}
	
	return block;
}

double *generateLinearFromTwoBlocks(int blocksize, int i, int j, int k, double **A, double **B) {
	double *block = calloc(blocksize * blocksize * 2, sizeof(double));
	
	int blockI, blockJ, blockK;
	int index = 0;
	for (blockI = i * blocksize; blockI < (i * blocksize) + blocksize; blockI++) {
		for (blockK = k * blocksize; blockK < (k * blocksize) + blocksize; blockK++) {
			block[index++] = A[blockI][blockK];	
		}
	}

	for (blockK = k * blocksize; blockK < (k * blocksize) + blocksize; blockK++) {
                for (blockJ = j * blocksize; blockJ < (j * blocksize) + blocksize; blockJ++) {
                        block[index++] = B[blockK][blockJ];     
                }
        }

	return block;
}


double **mult(int id, int blocksize, double *block, double **blockA, double **blockB, double **blockC) {
	int i, j, k;
	int index = 0;
	double sum;
	
	for (i = 0; i < blocksize; i++) {
		for (j = 0; j < blocksize; j++) {
			blockA[i][j] = block[index++];
		}
	}	

	for (i = 0; i < blocksize; i++) {
                for (j = 0; j < blocksize; j++) {
                        blockB[i][j] = block[index++];
                }
        }

	for (i = 0; i < blocksize; i++) {
		for (j = 0; j < blocksize; j++) {
			sum = 0.0;
			for (k = 0; k < blocksize; k++) {
				sum += blockA[i][k] * blockB[k][j];
			}
			blockC[i][j] = sum;
		}
	}
	return blockC;	
}


void usage(void) {
	fprintf(stderr, "usage: mmmw <blocksize> <N> <M> <P> <output file>\n");
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
