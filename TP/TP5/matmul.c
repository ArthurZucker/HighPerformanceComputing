/*
 * Sorbonne Universit�
 *
 * Programme de multiplication de matrices carr�es.
 */

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h>

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

#define REAL_T float
#define NB_TIMES 10

/*** Matmul: ***/
/* C += A x B
 * square matrices of order 'n'
 */
void matmul(int n, REAL_T * A, REAL_T * B, REAL_T * C)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}	/* for k */
		}		/* for j */
	}			/* for i */

}

void matmul_openmp(int n, REAL_T * A, REAL_T * B, REAL_T * C)
{
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				C[i * n + j] += A[i * n + k] * B[k * n + j];
			}	/* for k */
		}		/* for j */
	}			/* for i */

}

int main(int argc, char **argv)
{
	REAL_T *A, *B, *C;
	int n = 2;		/* default value */
	int nb = 0;

	/* Read 'n' on command line: */
	if (argc == 2) {
		n = atoi(argv[1]);
	}

	/* Allocate the matrices: */
	if ((A = (REAL_T *) malloc(n * n * sizeof(REAL_T))) == NULL)
		fprintf(stderr, "Error while allocating A.\n");
	if ((B = (REAL_T *) malloc(n * n * sizeof(REAL_T))) == NULL)
		fprintf(stderr, "Error while allocating B.\n");
	if ((C = (REAL_T *) malloc(n * n * sizeof(REAL_T))) == NULL)
		fprintf(stderr, "Error while allocating C.\n");

	/* Initialize the matrices */
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			*(A + i * n + j) = 1 / ((REAL_T) (i + j + 1));
			*(B + i * n + j) = 1.0;
			*(C + i * n + j) = 1.0;
		}

	/* Start timing */
	double debut = my_gettimeofday();
	for (nb = 0; nb < NB_TIMES; nb++) {
		/* Do matrix-product C=A*B+C */
		matmul(n, A, B, C);
		/* End timing */
	}
	double fin = my_gettimeofday();

	fprintf(stdout, "For n=%d: total computation time (with gettimeofday()) : %g s\n", n, (fin - debut) / NB_TIMES);
	fprintf(stdout, "For n=%d: performance = %g Gflop/s \n", n, (((double)2) * n * n * n / ((fin - debut) / NB_TIMES)) / ((double)1e9));	/* 2n^3 flops */

	/* Start timing */
	debut = my_gettimeofday();
	for (nb = 0; nb < NB_TIMES; nb++) {
		/* Do matrix-product C=A*B+C */
		matmul_openmp(n, A, B, C);
		/* End timing */
	}
	fin = my_gettimeofday();

	fprintf(stdout, "For n=%d: total computation time (with gettimeofday()) : %g s\n", n, (fin - debut) / NB_TIMES);
	fprintf(stdout, "For n=%d: performance = %g Gflop/s \n", n, (((double)2) * n * n * n / ((fin - debut) / NB_TIMES)) / ((double)1e9));	/* 2n^3 flops */

	/* Print 2x2 top-left square of C : */
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++)
			printf("%+e  ", C[i * n + j]);
		printf("\n");
	}
	printf("\n");

	/* Print 2x2 bottom-right square of C : */
	for (int i = n - 2; i < n; i++) {
		for (int j = n - 2; j < n; j++)
			printf("%+e  ", C[i * n + j]);
		printf("\n");
	}

	return EXIT_SUCCESS;
}
