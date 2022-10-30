#ifndef FVECTOR_H
#define FVECTOR_H


/* Structs definitions */

typedef struct fcsr {
	
	size_t n;  // Row count
	size_t m;  // Col count
	size_t *x; // Row indices
	size_t *y; // Col indices
	double *val;

} fcsr;

typedef struct fmatrix {
	
	size_t n; // Row count
	size_t m; // Col count
	double **val;
	
} fmatrix;

typedef struct fvector {

	size_t n;
	double *val;

} fvector;


/* Functions prototypes */

fcsr *init_fcsr(size_t n, size_t m);

fmatrix *init_fmatrix(size_t n, size_t m);

fvector *init_fvector(size_t n);

void free_csr(fcsr *my_csr);

void free_fmatrix(fmatrix *my_matrix);

void free_fvector(fvector *my_vector);

fvector *to_fvector(double *vals, size_t n);

double fvnorm(fvector *v, unsigned int p);

double fmnorm(fmatrix *mat, unsigned int p);

fvector *fones(size_t n);

fmatrix *feye(size_t n);

fmatrix *fdiag(fvector *v);

double fvdot(fvector *v1, fvector *v2);

#endif
