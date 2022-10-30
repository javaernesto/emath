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

	my_fcsr->x = NULL;
	my_fcsr->y = NULL;
	my_fcsr->val = NULL;
	
	free(my_fcsr);
}

/**
 * @brief - frees fvector structure
 * @param v - fvector pointer to free
 * @return void
 */

void free_fvector(fvector *v) {

	free(v->val);

	// for security reasons
	v->val = NULL;

	free(v);

	#ifdef FREE
    printf("FREE FV\n");
    #endif
}

/**
 * @brief - frees fmatrix structure
 * @param mat - fmatrix pointer to free
 * @return void
 */

void free_fmatrix(fmatrix *mat) {
	
	for (size_t i = 0; i < mat->n; i++)
	{
		free(mat->val[i]);

		// for security reasons
		mat->val[i] = NULL;
	}
	free(mat->val);

	// for security reasons
	mat->val = NULL;
	
	free(mat);

	#ifdef FREE
    printf("FREE FM\n");
    #endif
}

/**
 * @brief Auxiliary function for print_fvector
 * @param v integer vector to print
 */
void _print_fv(fvector *v) {
	
	printf("[");
	if (v->n <= 20) {
		for (size_t i = 0; i < v->n - 1; i++)
			printf("%lf, ", v->val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			printf("%lf, ", v->val[i]);
		printf( "..., ");

		for (size_t i = v->n - 3; i < v->n - 1; i++)
			printf("%lf, ", v->val[i]);
	}
	printf("%lf]", v->val[v->n - 1]);

}

/**
 * @brief Auxiliary function for fprint_fvector
 * @param v integer vector to print
 */
void _fprint_fv(FILE *fp, fvector *v) {
	
	fprintf(fp, "[");
	if (v->n <= 20) {
		for (size_t i = 0; i < v->n - 1; i++)
			fprintf(fp, "%lf, ", v->val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			fprintf(fp, "%lf, ", v->val[i]);
		fprintf(fp,  "..., ");

		for (size_t i = v->n - 3; i < v->n - 1; i++)
			fprintf(fp, "%lf, ", v->val[i]);
	}
	fprintf(fp, "%lf]", v->val[v->n - 1]);

}

/**
 * @brief Auxiliary function for print_fmatrix
 * @param val double int array containing row
 * @param n size of val
 */
void _print_lf(double *val, size_t n) {

	printf("[");
	if (n <= 21) {
		for (size_t i = 0; i < n - 1; i++)
			printf("%0.3lf, ", val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			printf("%0.3lf, ", val[i]);
		printf( "..., ");

		for (size_t i = n - 3; i < n - 1; i++)
			printf("%0.3lf, ", val[i]);
	}
	printf("%0.3lf]", val[n - 1]);
}

/**
 * @brief Auxiliary function for fprint_fmatrix
 * @param val double int array containing row
 * @param n size of val
 */
void _fprint_lf(FILE *fp, double *val, size_t n) {

	fprintf(fp, "[");
	if (n <= 21) {
		for (size_t i = 0; i < n - 1; i++)
			fprintf(fp, "%0.3lf, ", val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			fprintf(fp, "%0.3lf, ", val[i]);
		fprintf(fp,  "..., ");

		for (size_t i = n - 3; i < n - 1; i++)
			fprintf(fp, "%0.3lf, ", val[i]);
	}
	fprintf(fp, "%0.3lf]", val[n - 1]);
}

/**
 * @brief Print fvector (Python numpy library format). If vector size larger 
 * than 20, print brief representation
 * @param v double vector to print
 */
void print_fvector(fvector *v) {

	#ifdef DEBUG
	printf("size %ld\n", v->n);
	#endif

	printf("vector(");
	_print_fv(v);
	printf(")\n");
}

/**
 * @brief Print fvector to file (Python numpy library format). 
 * If vector size larger than 20, print brief representation
 * @param v double vector to print
 */
void fprint_fvector(FILE *fp, fvector *v) {

	#ifdef DEBUG
	fprintf(fp, "size %ld\n", v->n);
	#endif

	fprintf(fp, "vector(");
	_fprint_fv(fp, v);
	fprintf(fp, ")\n");
}

/**
 * @brief Print fmatrix (Python numpy library format). If row size or column 
 * size larger than 20, print brief representation
 * @param A double matrix to print
 */
void print_fmatrix(fmatrix *A) {

	#ifdef DEBUG
	printf("shape (%ld, %ld)\n", A->n, A->m);
	#endif

	printf("matrix([");
	if (A->n <= 21) {
		_print_lf(A->val[0], A->m);
		printf(",\n");	
		for (size_t i = 1; i < A->n - 1; i++) {
			printf("	");
			_print_lf(A->val[i], A->m);
			printf(",\n");
		}
		printf("	");
		_print_lf(A->val[A->n - 1], A->m);
	} else {
		_print_lf(A->val[0], A->m);
		printf(",\n");
		for (size_t i = 1; i < 3; i++) {
			printf("	");
			_print_lf(A->val[i], A->m);
			printf(",\n");
		}
		printf( "	...,\n");

		for (size_t i = A->n - 3; i < A->n - 1; i++) {
			printf("	");
			_print_lf(A->val[i], A->m);
			printf(",\n");
		}
		printf("	");
		_print_lf(A->val[A->n - 1], A->m);
	}
	printf("])\n");
}

/**
 * @brief Print fmatrix to file (Python numpy library format). 
 * If row size or column size larger than 20, print brief representation
 * @param A double matrix to print
 */
void fprint_fmatrix(FILE *fp, fmatrix *A) {

	#ifdef DEBUG
	fprintf(fp, "shape (%ld, %ld)\n", A->n, A->m);
	#endif

	fprintf(fp, "matrix([");
	if (A->n <= 21) {
		_fprint_lf(fp, A->val[0], A->m);
		fprintf(fp, ",\n");	
		for (size_t i = 1; i < A->n - 1; i++) {
			fprintf(fp, " 	");
			_fprint_lf(fp, A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp, "	");
		_fprint_lf(fp, A->val[A->n - 1], A->m);
	} else {
		_fprint_lf(fp, A->val[0], A->m);
		fprintf(fp, ",\n");
		for (size_t i = 1; i < 3; i++) {
			fprintf(fp, " 	");
			_fprint_lf(fp, A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp,  " 	...,\n");

		for (size_t i = A->n - 3; i < A->n - 1; i++) {
			fprintf(fp, "	");
			_fprint_lf(fp, A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp, " ");
		_fprint_lf(fp, A->val[A->n - 1], A->m);
	}
	fprintf(fp, "])\n");
}

/**
 * @brief Read csv or dat file (one column file only, use for structure.dat)
 * @param filename name of csv or dat file to read
 * @return fvector* 
 */
fmatrix *read_fcsv(char *filename) {

	FILE *fp;
	char line[1024];
	unsigned int lines = 0;
	size_t j = 0;
	int counter = 0;
	int ncol = 1;
    fmatrix *output;

	fp = fopen(filename, "r");
 
    if (!fp) {
        printf("Can't open file\n");
        
        return NULL;
    } else {
		
        // Check we are reading at least one column
        assert(ncol > 0);
		
		while (fgets(line, sizeof(line), fp)) {
			lines++;
		}

		/* Find number of separators to deduce number of columns */
		// TODO: Error handling. Inconsistent ncol
		while(line[j] != '\0') {
			if(line[j] == ',')
				counter++;
			j++;
		}

		ncol = counter + 1;
				
		/* Reset pointer to beginning of file after reading number of
		 * lines */
		fseek(fp, 0, SEEK_SET);

		output = init_fmatrix(lines, ncol);
				
		for (size_t i = 0; i < lines; i++) {
            for (size_t j = 0; j < ncol; j++)
			    fscanf(fp, "%lf%*[,\n]", &output->val[i][j]);
		}
		
		fclose(fp);
		
		return output;
	}
}

/**
 * @brief Print double vector to file
 * @param v double vector to save in file
 * @param filename name of file
 * @return int 0 if OK, 1 if cannot open file
 */
int fvtofile(fvector *v, const char *filename) {

	FILE *fp;
	size_t i = 0;

	fp = fopen(filename, "w");

	if (!fp) {
		printf("Cannot open file\n");

		return 1;
	} else {
		for (i = 0; i < v->n; i++)
			fprintf(fp, "%lf\n", v->val[i]);

		fclose(fp);

		return 0;
	}
}

/**
 * @brief Print double matrix to file
 * @param v double matrix to save in file
 * @param filename name of file
 * @return int 0 if OK, 1 if cannot open file
 */
int fmtofile(fmatrix *M, const char *filename) {

	FILE *fp;
	size_t i = 0, j = 0;

	fp = fopen(filename, "w");

	if (!fp) {
		printf("Cannot open file %s\n", filename);

		return 1;
	} else {
		for (i = 0; i < M->n; i++) {
			for (j = 0; j < M->m - 1; j++)
				fprintf(fp, "%lf,", M->val[i][j]);
			fprintf(fp, "%lf\n", M->val[i][j]);
		}

		fclose(fp);

		return 0;
	}
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
 * @brief Compute matrix of absolute values
 * @param A double matrix
 * @return fmatrix* 
 */
fmatrix *fmabs(fmatrix *A) {

	size_t i = 0, j = 0;
	fmatrix *output;

	output = init_fmatrix(A->n, A->m);

	for (i = 0; i < A->n; i++) {
		for (j = 0; j < A->m; j++) {
			output->val[i][j] = fabs(A->val[i][j]);
		}
	}

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
                sz += fabs((v->val)[i]);

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
                    sz += fabs(((mat->val)[i])[j]);
					
					if (sz > norm)
						norm = sz;
				}
				sz = 0;
			}
			break;
			
		case 1: // p = 1
			for (size_t j = 0; j < mat->n; j++) {
				for (size_t i = 0; i < mat->m; i++) {
                    sz += fabs(((mat->val)[i])[j]);
					
					if (sz > norm)
						norm = sz;
				}
				sz = 0;
			}
			break;

		case 2: // p = 2 (Frobenius)
			for (size_t j = 0; j < mat->n; j++) {
				for (size_t i = 0; i < mat->m; i++)
                    sz += mat->val[i][j] * mat->val[i][j];
			}
			norm = sqrt(sz);
			break;
			
		default:
			break;
	}
	
	return norm;
}

/**
 * @brief Compute p-distance between matrix A and matrix B
 * @param A double matrix
 * @param B double matrix
 * @param p 0 < p <= Inf such that we get ||A - B||_p (enter 0 for p = Inf)
 * @return double 
 */
double fmdist(fmatrix *A, fmatrix *B, unsigned int p) {

	double dist = 0.0;
	fmatrix *tmp, *C;

	/* Compute A - B first and then |A - B| */
	tmp = fmsub(A, B);
	C = fmabs(tmp);
	dist = fmnorm(C, p);

	free_fmatrix(tmp);
	free_fmatrix(C);

	return dist;
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
/**
 * @brief return matrix addition of two double matrices
 * @param A double matrix A
 * @param B double matrix B
 * @return *fmatrix
 */
fmatrix *fmadd(fmatrix *A, fmatrix *B) {

    // Assert fmatrix have same dimensions
	assert((A->n == B->n) && (A->m == B->m));

    fmatrix *C = init_fmatrix(A->n, A->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            (((C->val)[i]))[j] = (((A->val)[i]))[j] + (((B->val)[i]))[j];
        }
    }

    return C;
}

/**
 * @brief return matrix substraction of two double matrices
 * @param A double matrix A
 * @param B double matrix B
 * @return *fmatrix
 */
fmatrix *fmsub(fmatrix *A, fmatrix *B) {

    // Assert fmatrix have same dimensions
	assert((A->n == B->n) && (A->m == B->m));

    fmatrix *C = init_fmatrix(A->n, A->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            (((C->val)[i]))[j] = (((A->val)[i]))[j] - (((B->val)[i]))[j];
        }
    }

    return C;
}

/**
 * @brief return matrix division A / B (term by term) of two double matrices
 * @param A double matrix A
 * @param B double matrix B
 * @return *fmatrix
 */
fmatrix *fmdiv(fmatrix *A, fmatrix *B) {

    // Assert fmatrix have same dimensions
	assert((A->n == B->n) && (A->m == B->m));

    fmatrix *C = init_fmatrix(A->n, A->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            (((C->val)[i]))[j] = (((A->val)[i]))[j] / (((B->val)[i]))[j];
        }
    }

    return C;
}

/**
 * @brief return double matrix multiplication A * B
 * @param A double matrix
 * @param B double matrix
 * @return *fmatrix
 */
fmatrix *fmmul(fmatrix *A, fmatrix *B) {

    // Assert dimensions are compatible for multiplication
    assert((A->m == B->n));

    double c_ij = 0;
    fmatrix *C = init_fmatrix(A->n, B->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            for (size_t k = 0; k < A->m; k++) {
                c_ij += A->val[i][k] * B->val[k][j];
            }
            C->val[i][j] = c_ij;
            c_ij = 0;
        }
    }

    return C;
}

/**
 * @brief Compute coordinates of maximum element in double matrix A
 * @param A double matrix
 * @return size_t* 
 */
size_t *fmmax(fmatrix *A) {

	double counter = 0;
	
	size_t *coor = calloc(2, sizeof(size_t));

	for (int i = 0; i < A->n; i++) {
		for (int j = 0; j < A->m; j++) {
			if (A->val[i][j] > counter) {
				counter = A->val[i][j];
				coor[0] = i;
				coor[1] = j;
			}
		}
	}
	
	return coor;
}

/**
 * @brief Compute coordinates of minimum element in double matrix A
 * @param A double matrix
 * @return size_t* 
 */
size_t *fmmin(fmatrix *A) {

	double counter = 1;
	
	size_t *coor = calloc(2, sizeof(size_t));

	for (int i = 0; i < A->n; i++) {
		for (int j = 0; j < A->m; j++) {
			if (A->val[i][j] < counter) {
				counter = A->val[i][j];
				coor[0] = i;
				coor[1] = j;
			}
		}
	}
	
	return coor;
}

/**
 * @brief Find element candidate in vector v and return vector of same size 
 * containing 1 if candidate is in position i and 0 otherwise
 * @param v double vector to search
 * @param candidate element to search in vector
 * @return ivector* 
 */
ivector *fvfind(fvector *v, double candidate) {

	size_t k = 0;
	ivector *indices;

	indices = init_ivector(v->n);
	 
	for (k = 0; k < v->n; k++) {
		if (v->val[k] == candidate) {
			indices->val[k] = 1;		
		} else {
			indices->val[k] = 0;
		}
	}
	
	return indices;
}

/**
 * @brief Transform single-column double matrix to an double vector. Free 
 * matrix pointer.
 * @param A matrix to vectorize
 * @return *fvector
 */
fvector *fvectorize(fmatrix *A) {

	fvector *output;

	// Check matrix has one column only
	assert(A->m == 1);

	output = init_fvector(A->n);

	for (size_t i = 0; i < A->n; i++)
		output->val[i] = A->val[i][0];

	free_fmatrix(A);

	return output;
}


/**
 * @brief Transpose double matrix A
 * @param A double matrix to transpose
 * @return fmatrix*
 */
fmatrix *ftranspose(fmatrix *A) {

	fmatrix *B = init_fmatrix(A->m, A->n);

	for (size_t i = 0; i < A->n; i++) {
		for (size_t j = 0; j < A->m; j++) {
			B->val[i][j] = A->val[j][i];
		}
	}

	return B;
}