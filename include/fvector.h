#ifndef FVECTOR_H
#define FVECTOR_H

#include "ivector.h"


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

void free_fcsr(fcsr *my_csr);

void free_fmatrix(fmatrix *my_matrix);

void free_fvector(fvector *my_vector);

void _print_fv(fvector *v);

void _fprint_fv(FILE *fp, fvector *v);

void _print_lf(double *val, size_t n);

void _fprint_lf(FILE *fp, double *val, size_t n);

void print_fvector(fvector *v);

void fprint_fvector(FILE *fp, fvector *v);

void print_fmatrix(fmatrix *A);

void fprint_fmatrix(FILE *fp, fmatrix *A);

fmatrix *read_fcsv(char *filename);

int fvtofile(fvector *v, const char *filename);

int fmtofile(fmatrix *M, const char *filename);

fvector *to_fvector(double *vals, size_t n);

fmatrix *fmabs(fmatrix *A);

double fvnorm(fvector *v, unsigned int p);

double fmnorm(fmatrix *mat, unsigned int p);

double fmdist(fmatrix *A, fmatrix *B, unsigned int p);

fvector *fones(size_t n);

fmatrix *feye(size_t n);

fmatrix *fdiag(fvector *v);

double fvdot(fvector *v1, fvector *v2);

fmatrix *fmadd(fmatrix *A, fmatrix *B);

fmatrix *fmsub(fmatrix *A, fmatrix *B);

fmatrix *fmdiv(fmatrix *A, fmatrix *B);

fmatrix *fmmul(fmatrix *A, fmatrix *B);

size_t *fmmax(fmatrix *A);

size_t *fmmin(fmatrix *A);

ivector *fvfind(fvector *v, double candidate);

fvector *fvectorize(fmatrix *A);

fmatrix *ftranspose(fmatrix *A);

#endif
