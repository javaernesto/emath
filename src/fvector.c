#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "../include/fvector.h"


/**
 * @brief - initializes fcsr structure
 * @param n - number of rows
 * @param m - ajusted number of cols
 * @return *fcsr 
 */

fcsr *init_fcsr( size_t n, size_t m) {

	fcsr *my_fcsr = malloc(sizeof(struct fcsr));
	my_fcsr->n = n;
	my_fcsr->m = m;

	// Initialize values and indices
	my_fcsr->x = calloc(n, sizeof(size_t));
	my_fcsr->y = calloc(m, sizeof(size_t));

    my_fcsr->val = calloc(m, sizeof(double));	
	
	return my_fcsr;
}

/**
 * @brief - initializes fvector structure
 * @param n - size of fvector
 * @return *fvector
 */

fvector *init_fvector(size_t n) {

	fvector *my_fvector = malloc(sizeof(struct fvector));
	my_fvector->n = n;

	my_fvector->val = calloc(n, sizeof(double));	

	return my_fvector;
}

/**
 * @brief - initializes fmatrix structure
 * @param n - number of rows
 * @param m - number of columns
 * @return *fmatrix
 */

fmatrix *init_fmatrix(size_t n, size_t m) {
	
	/* Initialize ffmatrix structure to 0 */
	
	fmatrix *my_fmatrix = malloc(sizeof(struct fmatrix));
	my_fmatrix->n = n;
	my_fmatrix->m = m;

    my_fmatrix->val = malloc(n * sizeof(double *));
    
    for (size_t i = 0; i < n; i++) {
        my_fmatrix->val[i] = calloc(m, sizeof(double));
    }

	return my_fmatrix;
}

/**
 * @brief - frees fcsr structure
 * @param my_fcsr - fcsr pointer to free
 * @return void
 */

void free_fcsr(fcsr *my_fcsr) {
	
	/* Free fcsr structure */
	
	free(my_fcsr->x);
	free(my_fcsr->y);
	free(my_fcsr->val);
	
	free(my_fcsr);
}

/**
 * @brief - frees fvector structure
 * @param v - fvector pointer to free
 * @return void
 */

void free_fvector(fvector *v) {

	free(v->val);

	free(v);
}

/**
 * @brief - frees fmatrix structure
 * @param mat - fmatrix pointer to free
 * @return void
 */

void free_fmatrix(fmatrix *mat) {
	
	for (size_t i = 0; i < mat->n; i++)
		free(mat->val[i]);
	free(mat->val);
	
	free(mat);
}

/**
 * @brief - converts array vals of type t to fvector
 * @param vals - array of data
 * @param n - size of vals
 * @return *fvector 
 */

fvector *to_fvector(double *vals, size_t n) {
	
	fvector *output = init_fvector(n);

	memcpy(output->val, vals, n * sizeof(double));
	
	return output;
}

/**
 * @brief - return p-norm of fvector v
 * @param v - fvector v
 * @param p - 0 < p <= Inf such that we get ||v||_p (enter 0 for p = Inf)
 * @return double 
 */

double fvnorm(fvector *v, unsigned int p) {
	
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
 * @brief - return induced p-norm of fmatrix mat
 * @param mat - fmatrix mat
 * @param p - 0 < p <= Inf such that we get ||mat||_p (enter 0 for p = Inf)
 * @return double 
 */

double fmnorm(fmatrix *mat, unsigned int p) {

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
 * @brief - return fvector of ones of size n and type t
 * @param n - size n
 * @return *fvector
 */
fvector *fones(size_t n) {

	fvector *v = init_fvector(n);

	for (size_t i = 0; i < n; i++)
        (v->val)[i] = 1;

	return v;
}

/**
 * @brief - return identity fmatrix of size n 
 * @param n - size of identity fmatrix
 * @return *fmatrix 
 */
fmatrix *feye(size_t n) {

	fmatrix *mat = init_fmatrix(n, n);

    for (size_t i = 0; i < n; i++)
		(((mat->val)[i]))[i] = 1;

	return mat;
}

/**
 * @brief - return diagonal fmatrix with fvector v on diagonal
 * @param v - fvector to appear on the diagonal
 * @return *fmatrix
 */
fmatrix *fdiag(fvector *v ) {

	fmatrix *mat = init_fmatrix(v->n, v->n);

    for (size_t i = 0; i < mat->n; i++) {
        (((mat->val)[i]))[i] = (v->val)[i];
    }

	return mat;
}

/**
 * @brief - return dot product of two fvectors
 * @param a - fvector a
 * @param b - fvector b
 * @return double 
 */
double fvdot(fvector *v1, fvector *v2) {

	double prod = 0;

	// Assert fvector have same size
	assert(v1->n == v2->n);

	for (size_t i = 0; i < v1->n; i++)
        prod += (v1->val)[i] * (v2->val)[i];

	return prod;
}

int main() {

    // TESTING

	double v1[5] = {8.1, 2.4, 3.1, 4.9, 5.2};
	double v2[5] = {3, 0, 4.55, 0, 0};

	fvector *a = to_fvector(v1, 5);
	fvector *b = to_fvector(v2, 5);

	printf("ivdot %f \n", fvdot(a, b));

	fmatrix *mat = fdiag(a);
	for (size_t i = 0; i < mat->n; i++) {
		for (size_t j = 0; j < mat->m; j++)
			printf("%f ", (((mat->val)[i]))[j]);
		printf("\n");
	}

    printf("norm a = %.3f\n", fvnorm(a, 0));
    printf("norm b = %.3f\n", fvnorm(b, 2));

	free_fvector(a);
    free_fvector(b);
	free_fmatrix(mat);

	return 0;
}
