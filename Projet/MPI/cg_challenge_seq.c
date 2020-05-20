/* 
 * Sequential implementation of the Conjugate Gradient Method.
 *
 * Authors : Lilia Ziane Khodja & Charles Bouillaguet
 *
 * v1.01-challenge (2020-05-18)
 *
 * CHANGE LOG
 *  
 * USAGE: 
 * 	$ ./cg_challenge
 *
 * This code is almost identical to the "normal" cg.c.... 
 * EXCEPT THAT:
 *  + most integers are now 64 bits (cf. typedef ... i64 below)
 *  + the matrix is not loaded from a file
 *  + the matrix is pseudo-randomly generated build the build_mm() function()
 *  + the first argument of build_mm() is the size of the matrix. Can be adjusted.
 *  + sp_gemv and cg_solve are not modified (except for the 64-bit integers)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <inttypes.h>

#define THRESHOLD 1e-8		// maximum tolerance threshold

typedef int64_t i64;            // need 64-bit ints for more than 2^32 entries

struct csr_matrix_t {
	i64 n;			// dimension (64-bit)
	i64 nz;			// number of non-zero entries (64-bit)
	i64 *Ap;		// row pointers (64-bit)
	i64 *Aj;		// column indices (64-bit)
	double *Ax;		// actual coefficient
};

/*************************** Utility functions ********************************/

/* Seconds (wall-clock time) since an arbitrary point in the past */
double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

/* Pseudo-random function to initialize b (rumors says it comes from the NSA) */
#define ROR(x, r) ((x >> r) | (x << (64 - r)))
#define ROL(x, r) ((x << r) | (x >> (64 - r)))
#define R(x, y, k) (x = ROR(x, 8), x += y, x ^= k, y = ROL(y, 3), y ^= x)
double PRF(i64 i, unsigned long long seed)
{
	unsigned long long y = i, x = 0xBaadCafe, b = 0xDeadBeef, a = seed;
	R(x, y, b);
	for (i64 i = 0; i < 31; i++) {
		R(a, b, i);
		R(x, y, b);
	}
	x += i;
	union { double d; unsigned long long l;	} res;
	res.l = ((x << 2) >> 2) | (((1 << 10) - 1ll) << 52);
	return 2 * (res.d - 1.5);
}

/*************************** Matrix IO ****************************************/

/* generate a gaussian random variable. Pr(X == a) == exp(-a^2 / 2). Cf. wikipedia. */
static double normal_deviate(i64 i, i64 j)
{
	while (1) {
		double U = PRF(i, j);
		j = (2 * j + 13) % 0x7fffffff;
		double V = PRF(i, j);
		j = (2 * j + 13) % 0x7fffffff;
		double S = U*U + V*V;
		if (S >= 1)
			continue;
		return U * sqrt(-2 * log(S) / S);
	}
}

/* Generate a pseudo-random symetric positive definite matrix of size n. */
struct csr_matrix_t *build_mm(i64 n, double easyness)
{
	i64 nzmax = 64 * n; /* this is a crude upper-bound */

	/* -------- Directly assemble a CSR matrix ----- */
	double start = wtime();

	/* allocate CSR matrix */
	struct csr_matrix_t *A = malloc(sizeof(*A));
	if (A == NULL)
		err(1, "malloc failed");
	i64 *w = malloc((n + 1) * sizeof(*w));
	i64 *Ap = malloc((n + 1) * sizeof(*Ap));
	i64 *Aj = malloc(nzmax * sizeof(*Ap));
	double *Ax = malloc(nzmax * sizeof(*Ax));
	if (w == NULL || Ap == NULL || Aj == NULL || Ax == NULL)
		err(1, "Cannot allocate (CSR) sparse matrix");

	i64 k = 0;
	double scale_a = easyness * log2(n);
	double scale_b = log2(n);

	for (i64 i = 0; i < n; i++) {
		/* generate the i-th row of the matrix */

		/* diagonal entry */
		Ap[i] = k;
		Aj[k] = i;
		Ax[k] = scale_a + scale_b * normal_deviate(i, i);

		/* other entries of the i-th row (below and above diagonal) */
		i64 r = 1;
		for (i64 l = 1; l < n; l *= 2) {
			i64 u = i + l;
			i64 v = i - l;
			if (u < n) {
				Aj[k + r] = u;
				Ax[k + r] = -1 + normal_deviate(i*u, i+u);
				r++;
			}
			if (0 <= v) {
				Aj[k + r] = v;
				Ax[k + r] = -1 + normal_deviate(i*v, i+v);
				r++;
			}
		}
		k += r;
	}
	Ap[n] = k;

	double stop = wtime();
	fprintf(stderr, "     ---> nnz = %" PRId64 "\n", k);
	fprintf(stderr, "     ---> Assembled in CSR format in %.1fs\n", stop - start);
	fprintf(stderr, "     ---> CSR matrix size = %.1fMbyte\n", 1e-6 * (16. * k + 8. * n));

	A->n = n;
	A->nz = k;
	A->Ap = Ap;
	A->Aj = Aj;
	A->Ax = Ax;
	return A;
}

/*************************** Matrix accessors (unchanged) *********************************/

/* Copy the diagonal of A into the vector d. */
void extract_diagonal(const struct csr_matrix_t *A, double *d)
{
	i64 n = A->n;
	i64 *Ap = A->Ap;
	i64 *Aj = A->Aj;
	double *Ax = A->Ax;
	for (i64 i = 0; i < n; i++) {
		d[i] = 0.0;
		for (i64 u = Ap[i]; u < Ap[i + 1]; u++)
			if (i == Aj[u])
				d[i] += Ax[u];
	}
}

/* Matrix-vector product (with A in CSR format) : y = Ax */
void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y)
{
	i64 n = A->n;
	i64 *Ap = A->Ap;
	i64 *Aj = A->Aj;
	double *Ax = A->Ax;
	for (i64 i = 0; i < n; i++) {
		y[i] = 0;
		for (i64 u = Ap[i]; u < Ap[i + 1]; u++) {
			i64 j = Aj[u];
			double A_ij = Ax[u];
			y[i] += A_ij * x[j];
		}
	}
}

/*************************** Vector operations (unchaged) ********************************/

/* dot product */
double dot(const i64 n, const double *x, const double *y)
{
	double sum = 0.0;
	for (i64 i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

/* euclidean norm (a.k.a 2-norm) */
double norm(const i64 n, const double *x)
{
	return sqrt(dot(n, x, x));
}

/*********************** conjugate gradient algorithm (unchanged) *************************/

/* Solve Ax == b (the solution is written in x). Scratch must be preallocated of size 6n */
void cg_solve(const struct csr_matrix_t *A, const double *b, double *x, const double epsilon, double *scratch)
{
	i64 n = A->n;
	i64 nz = A->nz;

	fprintf(stderr, "[CG] Starting iterative solver\n");
	fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (16.0 * nz + 52.0 * n));
	fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n);

	double *r = scratch + n;	// residue
	double *z = scratch + 2 * n;	// preconditioned-residue
	double *p = scratch + 3 * n;	// search direction
	double *q = scratch + 4 * n;	// q == Ap
	double *d = scratch + 5 * n;	// diagonal entries of A (Jacobi preconditioning)

	/* Isolate diagonal */
	extract_diagonal(A, d);

	/* 
	 * This function follows closely the pseudo-code given in the (english)
	 * Wikipedia page "Conjugate gradient method". This is the version with
	 * preconditionning.
	 */

	/* We use x == 0 --- this avoids the first matrix-vector product. */
	for (i64 i = 0; i < n; i++)
		x[i] = 0.0;
	for (i64 i = 0; i < n; i++)	// r <-- b - Ax == b
		r[i] = b[i];
	for (i64 i = 0; i < n; i++)	// z <-- M^(-1).r
		z[i] = r[i] / d[i];
	for (i64 i = 0; i < n; i++)	// p <-- z
		p[i] = z[i];

	double rz = dot(n, r, z);
	double start = wtime();
	double last_display = start;
	int iter = 0;
	while (norm(n, r) > epsilon) {
		/* loop invariant : rz = dot(r, z) */
		double old_rz = rz;
		sp_gemv(A, p, q);	/* q <-- A.p */
		double alpha = old_rz / dot(n, p, q);
		for (i64 i = 0; i < n; i++)	// x <-- x + alpha*p
			x[i] += alpha * p[i];
		for (i64 i = 0; i < n; i++)	// r <-- r - alpha*q
			r[i] -= alpha * q[i];
		for (i64 i = 0; i < n; i++)	// z <-- M^(-1).r
			z[i] = r[i] / d[i];
		rz = dot(n, r, z);	// restore invariant
		double beta = rz / old_rz;
		for (i64 i = 0; i < n; i++)	// p <-- z + beta*p
			p[i] = z[i] + beta * p[i];
		iter++;
		double t = wtime();
		if (t - last_display > 0.5) {
			/* verbosity */
			double rate = iter / (t - start);	// iterations per s.
			double GFLOPs = 1e-9 * rate * (2 * nz + 12 * n);
			fprintf(stderr, "\r     ---> error : %2.2e, iter : %d (%.1f it/s, %.2f GFLOPs)", norm(n, r), iter, rate, GFLOPs);
			fflush(stdout);
			last_display = t;
		}
	}
	fprintf(stderr, "\n     ---> Finished in %.1fs and %d iterations\n", wtime() - start, iter);
}

/******************************* main program *********************************/

/* options descriptor */
struct option longopts[4] = {
	{"seed", required_argument, NULL, 's'},
	{"solution", required_argument, NULL, 'o'},
	{NULL, 0, NULL, 0}
};

int main(int argc, char **argv)
{
	/* Parse command-line options */
	long long seed = 0;
	char *solution_filename = NULL;
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 's':
			seed = atoll(optarg);
			break;
		case 'o':
			solution_filename = optarg;
			break;
		default:
			errx(1, "Unknown option");
		}
	}

	/* Build the matrix --- WARNING, THIS ALLOCATES 400GB! */
	struct csr_matrix_t *A = build_mm(450000000, 5);

	/* Allocate memory */
	i64 n = A->n;
	double *mem = malloc(8 * n * sizeof(double)); /* WARNING, THIS ALLOCATES 26GB. */
	if (mem == NULL)
		err(1, "cannot allocate dense vectors");
	double *x = mem;	/* solution vector */
	double *b = mem + n;	/* right-hand side */
	double *scratch = mem + 2 * n;	/* workspace for cg_solve() */

	/* Prepare right-hand size */
	for (i64 i = 0; i < n; i++)
		b[i] = PRF(i, seed);
	
	/* solve Ax == b */
	cg_solve(A, b, x, THRESHOLD, scratch);

	/* Dump the solution vector */
	FILE *f_x = stdout;
	if (solution_filename != NULL) {
		f_x = fopen(solution_filename, "w");
		if (f_x == NULL)
			err(1, "cannot open solution file %s", solution_filename);
		fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
	}
	for (i64 i = 0; i < n; i++)
		fprintf(f_x, "%a\n", x[i]);
	return EXIT_SUCCESS;
}