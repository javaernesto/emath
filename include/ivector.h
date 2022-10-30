#ifndef IVECTOR_H
#define IVECTOR_H


/* Structs definitions */

/**
 * @brief Struct for Column Sparse Representation (CSR) integer matrix
 */
typedef struct icsr {
	
	size_t n;  // Row count
	size_t m;  // Col count
	size_t *x; // Row indices
	size_t *y; // Col indices
	long *val;

} icsr;

typedef struct imatrix {
	
	size_t n; // Row count
	size_t m; // Col count
	long **val;
	
} imatrix;

typedef struct ivector {

	size_t n;
	long *val;

} ivector;


/* Functions prototypes */

icsr *init_icsr(size_t n, size_t m);

imatrix *init_imatrix(size_t n, size_t m);

ivector *init_ivector(size_t n);

void free_csr(icsr *my_csr);

void free_imatrix(imatrix *my_matrix);

void free_ivector(ivector *my_vector);

ivector *to_ivector(long *vals, size_t n);

double ivnorm(ivector *v, unsigned int p);

double imnorm(imatrix *mat, unsigned int p);

ivector *iones(size_t n);

imatrix *ieye(size_t n);

imatrix *idiag(ivector *v);

ivector *ivadd(ivector *a, ivector *b);

imatrix *imadd(imatrix *A, imatrix *B);

long ivdot(ivector *v1, ivector *v2);

long itrace(imatrix *A);

imatrix *immul(imatrix *A, imatrix *B);

#endif
