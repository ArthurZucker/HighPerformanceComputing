#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
//MPI => Utiliser des coeurs d'autes PC
//OpenMP => Utiliser tout les coeurs de son PC

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int** produit_matriciel(int n){
  int** C = (int**)calloc(n,sizeof(int*));
  int** A = (int**)calloc(n,sizeof(int*));
  int** B = (int**)calloc(n,sizeof(int*));
  for(int i=0;i<n;i++){
    C[i]=(int*)calloc(n,sizeof(int));
    A[i]=(int*)calloc(n,sizeof(int));
    B[i]=(int*)calloc(n,sizeof(int));
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<n;i++){
      A[i][j]= rand()%50;
      B[i][j]= rand()%50;
    }
  }

  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<n;k++){
        C[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
  return C;
}

int** produit_matriciel_openmp(int n){
  int** C = (int**)calloc(n,sizeof(int*));
  int** A = (int**)calloc(n,sizeof(int*));
  int** B = (int**)calloc(n,sizeof(int*));
  for(int i=0;i<n;i++){
    C[i]=(int*)calloc(n,sizeof(int));
    A[i]=(int*)calloc(n,sizeof(int));
    B[i]=(int*)calloc(n,sizeof(int));
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<n;i++){
      A[i][j]= rand()%50;
      B[i][j]= rand()%50;
    }
  }
  #pragma omp parallel for //mieux ici car on crée qu'une fois l'équipe de thread
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<n;k++){
        C[i][j]+=A[i][k]*B[k][j];
      }
    }
  }
  return C;
}

int main(int argc, char const *argv[]) {
  srand(time(NULL));
  int N=6;
  /* Chronometrage */
  double debut, fin;

  /* debut du chronometrage */
  debut = my_gettimeofday();

  produit_matriciel(N);

  /* fin du chronometrage */
  fin = my_gettimeofday();
  fprintf(stderr, "Temps total de calcul : %g sec\n", fin - debut);
  fprintf(stdout, "%g\n", fin - debut);

  /* debut du chronometrage */
  debut = my_gettimeofday();

  produit_matriciel_openmp(N);

  /* fin du chronometrage */
  fin = my_gettimeofday();
  fprintf(stderr, "Temps total de calcul : %g sec\n", fin - debut);
  fprintf(stdout, "%g\n", fin - debut);

  return 0;
}
