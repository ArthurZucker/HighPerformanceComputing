#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
//MPI => Utiliser des coeurs d'autes PC
//OpenMP => Utiliser tout les coeurs de son PC

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

double pi_seq(){
  double res = 0.0;
  long int N = 1000000000;

  for (int i=0;i<N;i++)
    res+=4/(1+(i/N)^2);

  double pi=res/N;
  return pi;
}

double pi_openmp(){
  double res = 0.0;
  long int N = 1000000000;

  #pragma omp parallel for reduction(+:res)
  for (int i=0;i<N;i++)
    res+=4/(1+(i/N)^2);

  double pi=res/N;
  return pi;
}

int main(int argc, char const *argv[]) {
  /* Chronometrage */
  double debut, fin;

  /* debut du chronometrage */
  debut = my_gettimeofday();

  double pi = pi_seq();
  printf("pi: %f\n",pi);

  /* fin du chronometrage */
  fin = my_gettimeofday();
  fprintf(stderr, "Temps total de calcul : %g sec\n", fin - debut);
  fprintf(stdout, "%g\n", fin - debut);

  /* debut du chronometrage */
  debut = my_gettimeofday();

  pi = pi_openmp();
  printf("pi_openmp: %f\n",pi);

  /* fin du chronometrage */
  fin = my_gettimeofday();
  fprintf(stderr, "Temps total de calcul : %g sec\n", fin - debut);
  fprintf(stdout, "%g\n", fin - debut);

  return 0;
}
