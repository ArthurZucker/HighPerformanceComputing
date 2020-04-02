#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

///////////////////////////////////////////////////////////////////////
/// D'apres http://openclassrooms.com/courses/le-tri-rapide-qsort : 
void echanger(int tableau[], int a, int b)
{
	int temp = tableau[a];
	tableau[a] = tableau[b];
	tableau[b] = temp;
}

void QuickSort(int tableau[], int debut, int fin)
{
	int gauche = debut - 1;
	int droite = fin + 1;
	
	/* Si le tableau est de longueur nulle, il n'y a rien à faire. */
	if (debut >= fin)
		return;

	const int pivot = tableau[debut];	// premier element choisi comme pivot 

	/* Sinon, on parcourt le tableau, une fois de droite à gauche, et une
	   autre de gauche à droite, à la recherche d'éléments mal placés,
	   que l'on permute. Si les deux parcours se croisent, on arrête. */
	while (1) {
		do
			droite--;
		while (tableau[droite] > pivot);
		do
			gauche++;
		while (gauche <= fin && tableau[gauche] <= pivot);

		if (gauche < droite)
			echanger(tableau, gauche, droite);
		else
			break;
	}

	/* On met le pivot à sa place: */
	echanger(tableau, debut, droite);

	/* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
	   supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
	   pour cela... la méthode quickSort elle-même ! */
	QuickSort(tableau, debut, droite - 1);
	QuickSort(tableau, droite + 1, fin);
}

///////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	int n = 27;		/* valeur par defaut pour la taille du tableau en : 2^n */

	/* Lecture de 'n' sur la ligne de commande : */
	if (argc == 2)
		n = atoi(argv[1]);

	/* Allocation memoire et initialisation: */
	int taille = (int)1 << n;
	int *tableau = (int *)malloc(taille * sizeof(int));
	srand(time(NULL));
	for (int i = 0; i < taille; i++)
		tableau[i] = rand();

	/* Declenchement chrono */
	double debut = my_gettimeofday();

	/* Calcul :  */
	QuickSort(tableau, 0, taille - 1);

	/* Arret chrono */
	double fin = my_gettimeofday();

	printf("Pour n=%d (2^n=%d): temps de calcul (avec gettimeofday()) : %g s\n", n, taille, fin - debut);

	/* Verification : */
	printf("Verification ... ");
	fflush(stdout);
	for (int i = 0; i < taille - 1; i++) {
		if (tableau[i] > tableau[i + 1]) {
			printf("### Erreur pour tableau[%d] = %d et tableau[%d] = %d \n", i, tableau[i], i + 1, tableau[i + 1]);
			exit(EXIT_FAILURE);
		}
	}
	printf(" ok.\n");

	return 0;
}
