#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

void usage(void);
double cpu_time(void);

int main(int argc, char *argv[]) {
	if (argc < 3 || argc > 4) usage();

	int my_id;
	int nproc;
	int msgsize = atoi(argv[1]);
	int trials = atoi(argv[2]);
	int bidirectional = (argc == 4) ? 1 : 0;;
	double ctime;
	double ctime1;
	double ctime2;
	double *send0 = (double *) malloc(msgsize * sizeof(double));
	double *recv0 = (double *) calloc(msgsize, sizeof(double));
	double *send1 = (double *) calloc(msgsize, sizeof(double));
	double *recv1 = (double *) calloc(msgsize, sizeof(double));
	MPI_Status status[4];
	MPI_Request request[4];


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);	


	// Store current cpu clock
	if (my_id == 0) {
		int i;
		for (i = 0; i < msgsize; i++) {
			send0[i] = drand48();
		}
		ctime1 = cpu_time();
        } else {
		int i;
		for (i = 0; i < msgsize; i++) {
			send1[i] = drand48();
		}
	}


	// Send messages
	int i;
	for (i = 0; i < trials; i++) {
		if (bidirectional == 0) {
			if (my_id == 0) {
				MPI_Isend(send0, msgsize, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request[0]);
				MPI_Wait(&request[0], &status[0]);
				MPI_Irecv(recv0, msgsize, MPI_DOUBLE, my_id + 1, 1, MPI_COMM_WORLD, &request[1]);
				MPI_Wait(&request[1], &status[1]);
			} else {
				MPI_Irecv(recv1, msgsize, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request[0]);
				MPI_Wait(&request[0], &status[0]);
				MPI_Isend(recv1, msgsize, MPI_DOUBLE, my_id - 1, 1, MPI_COMM_WORLD, &request[1]);
				MPI_Wait(&request[1], &status[1]);
			}

			
			if (my_id == 0) {
				int j;
				for (j = 0; j < msgsize; j++) {
					if (send0[j] != recv0[j]) {
						printf("ERROR 0: %d\n", j);
					}
				}
			}
			
		} else {
			if (my_id == 0) {
				// Send proc 0 data
				MPI_Isend(send0, msgsize, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request[0]);
				MPI_Wait(&request[0], &status[0]);

				// Recv proc 1 data
				MPI_Irecv(recv1, msgsize, MPI_DOUBLE, my_id + 1, 1, MPI_COMM_WORLD, &request[1]);
				MPI_Wait(&request[1], &status[1]);				

				// Send proc 1 data
				MPI_Isend(recv1, msgsize, MPI_DOUBLE, my_id + 1, 3, MPI_COMM_WORLD, &request[2]);
				MPI_Wait(&request[2], &status[2]);

				// Recv proc 0 data
				MPI_Irecv(recv0, msgsize, MPI_DOUBLE, my_id + 1, 2, MPI_COMM_WORLD, &request[3]);
				MPI_Wait(&request[3], &status[3]);

			} else {
				// Recv proc 0 data
				MPI_Irecv(recv0, msgsize, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request[0]);
                                MPI_Wait(&request[0], &status[0]);

				// Send proc 1 data
				MPI_Isend(send1, msgsize, MPI_DOUBLE, my_id - 1, 1, MPI_COMM_WORLD, &request[1]);
				MPI_Wait(&request[1], &status[1]);

				// Recv proc 1 data
				MPI_Irecv(recv1, msgsize, MPI_DOUBLE, my_id - 1, 3, MPI_COMM_WORLD, &request[3]);
                                MPI_Wait(&request[3], &status[3]);

				// Send proc 0 data
				MPI_Isend(recv0, msgsize, MPI_DOUBLE, my_id - 1, 2, MPI_COMM_WORLD, &request[2]);
				MPI_Wait(&request[2], &status[2]);				
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			int j;
			for (j = 0; j < msgsize; j++) {
				if (my_id == 0) {
					if (send0[j] != recv0[j]) {
						printf("ERROR BIDIRECTIONAL 0: %d\n", msgsize);
					}
				} else {
					if (send1[j] != recv1[j]) {
						printf("ERROR BIDIRECTIONAL 1: %d\n", msgsize);
					}
				}
			}
		}
	}


	// Process 0 print results to sdout 
	if (my_id == 0) {
                ctime2 = cpu_time();
        	ctime = (ctime2 - ctime1) / (double) trials;
		double rate;

		if (bidirectional == 0) {
			rate = (ctime == 0) ? msgsize : (double) msgsize / ctime * 2;
			printf("%d bytes (%d trials) %f time %f Bytes/s\n", msgsize, trials, ctime, rate);
        	} else {
			rate = (ctime == 0) ? msgsize : (double) msgsize / ctime * 4;
			printf("%d bytes (%d trials) %f time %f Bytes/s (bidirectional)\n", msgsize, trials, ctime, rate);
		}
	}
	
	MPI_Finalize();
	return 0;
}

/** usage
 * Print usage message to stderr
*/
void usage(void) {
	fprintf(stderr, "usage: msgbench <msgsize> <trials> [bidirectional\n");
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
 *John Burkardt
 * Parameters:
 * Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
double cpu_time (void) {
        double value;
        value = (double) clock () / (double) CLOCKS_PER_SEC;
        return value;
}
