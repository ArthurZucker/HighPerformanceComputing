#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <sys/time.h>

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

///////////////////////////////////////////////////////////////////////
/// WARNING: naive algorithm with worst operation count!
int fib2(int n)
{
	if (n < 2){
		return n;
	}
	else{
		int i, j;
		i = fib2(n - 1);
		j = fib2(n - 2);
		return i + j;
	}
}

int fib(int n)
{
	if (n < 2){
		return n;
	}
	else if(n<20){
		return fib2(n);
	}
	else {
		int i, j;
		/*Utiliser les threads que quand c'est utile*/
		#pragma omp task shared(i)
		i = fib(n - 1);
		#pragma omp task shared(j)
		j = fib(n - 2);
		#pragma omp taskwait
		return i + j;
	}
}

///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	int n = 45;		/* default value -> roughly 10 seconds of computation */

	/* Read 'n' on command line: */
	if (argc == 2)
		n = atoi(argv[1]);

	/* Start timing */
	double debut = my_gettimeofday();

	/* Do computation:  */
	int res;
	#pragma omp parallel
	#pragma omp single
	res = fib(n);

	/* End timing */
	double fin = my_gettimeofday();

	printf("fib(%d)=%d\n", n, res);
	printf("For n=%d: total computation time (with gettimeofday()) : %g s\n", n, fin - debut);

	return 0;
}
