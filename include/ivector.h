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

void free_icsr(icsr *my_csr);

void free_imatrix(imatrix *my_matrix);

void free_ivector(ivector *my_vector);

void _print_iv(ivector *v);

void _fprint_iv(FILE *fp, ivector *v);

void _print_ld(long *val, size_t n);

void _fprint_ld(FILE *fp, long *val, size_t n);

void print_ivector(ivector *v);

void fprint_ivector(FILE *fp, ivector *v);

void print_imatrix(imatrix *A);

void fprint_imatrix(FILE *fp, imatrix *A);

imatrix *read_icsv(char *filename);

int ivtofile(ivector *v, const char *filename);

int imtofile(imatrix *M, const char *filename);

ivector *to_ivector(long *vals, size_t n);

double ivnorm(ivector *v, unsigned int p);

double imnorm(imatrix *mat, unsigned int p);

ivector *iarange(size_t n);

ivector *iones(size_t n);

imatrix *ieye(size_t n);

imatrix *idiag(ivector *v);

ivector *ivadd(ivector *a, ivector *b);

ivector *ivsub(ivector *a, ivector *b);

imatrix *imadd(imatrix *A, imatrix *B);

imatrix *imsub(imatrix *A, imatrix *B);

size_t *immax(imatrix *A);

ivector *ivfind(ivector *v, long candidate);

long ivdot(ivector *v1, ivector *v2);

long itrace(imatrix *A);

imatrix *immul(imatrix *A, imatrix *B);

long int isum(ivector *v, size_t a, size_t b);

imatrix *itranspose(imatrix *A);

ivector *ivectorize(imatrix *A);

#endif
