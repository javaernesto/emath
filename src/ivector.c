#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ivector.h"


/**
 * @brief initialize icsr structure
 * @param n number of rows
 * @param m ajusted number of cols
 * @return *icsr 
 */
icsr *init_icsr( size_t n, size_t m) {

	icsr *my_icsr = malloc(sizeof(struct icsr));
	my_icsr->n = n;
	my_icsr->m = m;

	// Initialize values and indices
	my_icsr->x = calloc(n, sizeof(size_t));
	my_icsr->y = calloc(m, sizeof(size_t));

    my_icsr->val = calloc(m, sizeof(long));	
	
	return my_icsr;
}

/**
 * @brief initialize ivector structure
 * @param n size of ivector
 * @return *ivector
 */
ivector *init_ivector(size_t n) {

	ivector *my_ivector = malloc(sizeof(struct ivector));
	my_ivector->n = n;

	my_ivector->val = calloc(n, sizeof(long));	

	return my_ivector;
}

/**
 * @brief initialize imatrix structure
 * @param n number of rows
 * @param m number of columns
 * @return *imatrix
 */
imatrix *init_imatrix(size_t n, size_t m) {
	
	/* Initialize fimatrix structure to 0 */
	
	imatrix *my_imatrix = malloc(sizeof(struct imatrix));
	my_imatrix->n = n;
	my_imatrix->m = m;

    my_imatrix->val = malloc(n * sizeof(long *));
    
    for (size_t i = 0; i < n; i++) {
        my_imatrix->val[i] = calloc(m, sizeof(long));
    }

	return my_imatrix;
}

/**
 * @brief free icsr structure
 * @param my_icsr icsr pointer to free
 * @return void
 */
void free_icsr(icsr *my_icsr) {
	
	/* Free icsr structure */
	
	free(my_icsr->x);
	free(my_icsr->y);
	free(my_icsr->val);
	
	free(my_icsr);
}

/**
 * @brief free ivector structure
 * @param v ivector pointer to free
 * @return void
 */
void free_ivector(ivector *v) {

	free(v->val);

	free(v);
}

/**
 * @brief free imatrix structure
 * @param mat imatrix pointer to free
 * @return void
 */
void free_imatrix(imatrix *mat) {
	
	for (size_t i = 0; i < mat->n; i++)
		free(mat->val[i]);
	free(mat->val);
	
	free(mat);
}

/**
 * @brief convert integer array vals to ineteger vector
 * @param vals array of data
 * @param n size of vals
 * @return *ivector 
 */
ivector *to_ivector(long *vals, size_t n) {
	
	ivector *output = init_ivector(n);

	memcpy(output->val, vals, n * sizeof(long));
	
	return output;
}

/**
 * @brief return p-norm of integer vector v
 * @param v integer vector
 * @param p 0 < p <= Inf such that we get ||v||_p (enter 0 for p = Inf)
 * @return double 
 */
double ivnorm(ivector *v, unsigned int p) {
	
	double sz 	= 0.0;
	double norm = 0.0;

	switch(p) {
		case 0: // p = Inf -> max
			for (size_t i = 0; i < v->n; i++) {
				if ((v->val)[i] > norm) {
					norm = (v->val)[i];
                }
			}
			break;
			
		case 1: // p = 1	
			for (size_t i = 0; i< v->n; i++)
                sz += labs((v->val)[i]);

			norm = sz;
			break;

		case 2: // p = 2
			for (size_t i = 0; i < v->n; i++) 
            	sz += (v->val)[i] * (v->val)[i];

			norm = sqrt(sz);
			break;
			
		default:
			break;
	}
	
	return norm;
}

/**
 * @brief return induced p-norm of imatrix mat
 * @param mat integer matrix
 * @param p 0 < p <= Inf such that we get ||mat||_p (enter 0 for p = Inf)
 * @return double 
 */
double imnorm(imatrix *mat, unsigned int p) {

	double sz 	= 0.0;
	double norm = 0.0;

	switch(p) {
		case 0: // p = Inf -> max rows
			for (size_t i = 0; i < mat->n; i++) {
				for (size_t j = 0; j < mat->m; j++) {
                    sz += abs(((mat->val)[i])[j]);
					
					if (sz > norm)
						norm = sz;
				}
				sz = 0;
			}
			break;
			
		case 1: // p = 1
			for (size_t j = 0; j < mat->n; j++) {
				for (size_t i = 0; i < mat->m; i++) {
                    sz += abs(((mat->val)[i])[j]);
					
					if (sz > norm)
						norm = sz;
				}
				sz = 0;
			}
			break;
			
		default:
			break;
	}
	
	return norm;
}

/**
 * @brief return integer vector of ones of size 
 * @param n size of vector 
 * @return *ivector
 */
ivector *ones(size_t n) {

	ivector *v = init_ivector(n);

	for (size_t i = 0; i < n; i++)
        (v->val)[i] = 1;

	return v;
}

/**
 * @brief return identity integer matrix of size n 
 * @param n size of identity imatrix
 * @return *imatrix 
 */
imatrix *eye(size_t n) {

	imatrix *mat = init_imatrix(n, n);

    for (size_t i = 0; i < n; i++)
		(((mat->val)[i]))[i] = 1;

	return mat;
}

/**
 * @brief return diagonal integer matrix with ineteger vector v on diagonal
 * @param v integer vector to appear on the diagonal
 * @return *imatrix
 */
imatrix *diag(ivector *v) {

	imatrix *mat = init_imatrix(v->n, v->n);

    for (size_t i = 0; i < mat->n; i++) {
        (((mat->val)[i]))[i] = (v->val)[i];
    }

	return mat;
}

/**
 * @brief return dot product of two integer vectors
 * @param a integer vector a
 * @param b integer vector b
 * @return double 
 */
long vdot(ivector *v1, ivector *v2) {

	double prod = 0;

	// Assert ivector have same size
	assert(v1->n == v2->n);

	for (size_t i = 0; i < v1->n; i++)
        prod += (v1->val)[i] * (v2->val)[i];

	return prod;
}

/**
 * @brief return vector addition of two vectors
 * @param a integer vector a
 * @param b integer vector b
 * @return *ivector
 */
ivector *ivadd(ivector *a, ivector *b) {

    // Assert ivector have same size
	assert(a->n == b->n);

    ivector *c = init_ivector(a->n);

    for (size_t i = 0; i < c->n; i++)
        (c->val)[i] = (a->val)[i] + (b->val)[i];

    return c;
}

/**
 * @brief return matrix addition of two integer matrices
 * @param A integer matrix A
 * @param B integer matrix B
 * @return *imatrix
 */
imatrix *imadd(imatrix *A, imatrix *B) {

    // Assert imatrix have same dimensions
	assert((A->n == B->n) && (A->m == B->m));

    imatrix *C = init_imatrix(A->n, A->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            (((C->val)[i]))[j] = (((A->val)[i]))[j] + (((B->val)[i]))[j];
        }
    }

    return C;
}

/**
 * @brief return trace of square integer matrix
 * @param A integer square matrix 
 * @return long 
 */
long itrace(imatrix *A) {

    // Assert A matrix is square matrix
    assert(A->n == A->m);

    long trace = 0;

    for (size_t i = 0; i < A->n; i++) {
        trace += (((A->val)[i]))[i];
    }

    return trace;
}

/**
 * @brief return integer matrix multiplication A * B
 * @param A integer matrix
 * @param B integer matrix
 * @return *imatrix
 */
imatrix *immul(imatrix *A, imatrix *B) {

    // Assert dimensions are compatible for multiplication
    assert((A->m == B->n));

    long c_ij = 0;
    imatrix *C = init_imatrix(A->n, B->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            for (size_t k = 0; k < A->m; k++) {
                c_ij += (((A->val)[i]))[k] * (((B->val)[k]))[j];
            }
            (((C->val)[i]))[i] = c_ij;
            c_ij = 0;
        }
    }

    return C;
}

int main() {

    /* TESTING */

	long v1[5] = {8, 2, 3, 4, 5};
	long v2[5] = {3, 0, 4, 0, 0};

	ivector *a = to_ivector(v1, 5);
	ivector *b = to_ivector(v2, 5);

    ivector *c = ivadd(a, b);

	printf("vdot %ld \n", vdot(a, b));

    for (size_t i = 0; i < c->n; i++) {
        printf("%ld ", (c->val)[i]);
    }
    printf("\n");

	imatrix *mat = diag(a);
	for (size_t i = 0; i < mat->n; i++) {
		for (size_t j = 0; j < mat->m; j++)
			printf("%ld ", (((mat->val)[i]))[j]);
		printf("\n");
	}

    printf("norm a = %.3f\n", ivnorm(a, 0));
    printf("norm b = %.3f\n", ivnorm(b, 2));

	free_ivector(a);
    free_ivector(b);
    free_ivector(c);
	free_imatrix(mat);

	return 0;
}
