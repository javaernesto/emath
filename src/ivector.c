#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "../include/ivector.h"


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

	#ifdef FREE
    printf("FREE IV\n");
    #endif
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

	#ifdef FREE
    printf("FREE IM\n");
    #endif
}

/**
 * @brief Auxiliary function for print_ivector
 * @param v integer vector to print
 */
void _print_iv(ivector *v) {
	
	printf("[");
	if (v->n <= 20) {
		for (size_t i = 0; i < v->n - 1; i++)
			printf("%ld, ", v->val[i]);
	} else {
		for (size_t i = 0; i < 20; i++)
			printf("%ld, ", v->val[i]);
		printf( "..., ");

		for (size_t i = v->n - 20; i < v->n - 1; i++)
			printf("%ld, ", v->val[i]);
	}
	printf("%ld]", v->val[v->n - 1]);

}

/**
 * @brief Auxiliary function for fprint_ivector
 * @param v integer vector to print
 */
void _fprint_iv(FILE *fp, ivector *v) {
	
	fprintf(fp, "[");
	if (v->n <= 20) {
		for (size_t i = 0; i < v->n - 1; i++)
			fprintf(fp, "%ld, ", v->val[i]);
	} else {
		for (size_t i = 0; i < 20; i++)
			fprintf(fp, "%ld, ", v->val[i]);
		fprintf(fp,  "..., ");

		for (size_t i = v->n - 20; i < v->n - 1; i++)
			fprintf(fp, "%ld, ", v->val[i]);
	}
	fprintf(fp, "%ld]", v->val[v->n - 1]);

}

/**
 * @brief Auxiliary function for print_imatrix
 * @param val long int array containing row
 * @param n size of val
 */
void _print_ld(long *val, size_t n) {

	printf("[");
	if (n <= 20) {
		for (size_t i = 0; i < n - 1; i++)
			printf("%ld, ", val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			printf("%ld, ", val[i]);
		printf( "..., ");

		for (size_t i = n - 3; i < n - 1; i++)
			printf("%ld, ", val[i]);
	}
	printf("%ld]", val[n - 1]);
}

/**
 * @brief Auxiliary function for fprint_imatrix
 * @param val long int array containing row
 * @param n size of val
 */
void _fprint_ld(FILE *fp, long *val, size_t n) {

	fprintf(fp, "[");
	if (n <= 20) {
		for (size_t i = 0; i < n - 1; i++)
			fprintf(fp, "%ld, ", val[i]);
	} else {
		for (size_t i = 0; i < 3; i++)
			fprintf(fp, "%ld, ", val[i]);
		fprintf(fp,  "..., ");

		for (size_t i = n - 3; i < n - 1; i++)
			fprintf(fp, "%ld, ", val[i]);
	}
	fprintf(fp, "%ld]", val[n - 1]);
}

/**
 * @brief Print ivector (Python numpy library format). If vector size larger 
 * than 20, print brief representation
 * @param v integer vector to print
 */
void print_ivector(ivector *v) {

	#ifdef DEBUG
	printf("size %ld\n", v->n);
	#endif

	printf("vector(");
	_print_iv(v);
	printf(")\n");
}

/**
 * @brief Print ivector to file (Python numpy library format). 
 * If vector size larger than 20, print brief representation
 * @param v integer vector to print
 */
void fprint_ivector(FILE *fp, ivector *v) {

	#ifdef DEBUG
	fprintf(fp, "size %ld\n", v->n);
	#endif

	fprintf(fp, "vector(");
	_fprint_iv(fp, v);
	fprintf(fp, ")\n");
}

/**
 * @brief Print imatrix (Python numpy library format). If row size or column 
 * size larger than 20, print brief representation
 * @param A integer matrix to print
 */
void print_imatrix(imatrix *A) {

	#ifdef DEBUG
	printf("shape (%ld, %ld)\n", A->n, A->m);
	#endif

	printf("matrix([");
	if (A->n <= 20) {
		_print_ld(A->val[0], A->m);
		printf(",\n");	
		for (size_t i = 1; i < A->n - 1; i++) {
			printf("	");
			_print_ld(A->val[i], A->m);
			printf(",\n");
		}
		printf("	");
		_print_ld(A->val[A->n - 1], A->m);
	} else {
		_print_ld(A->val[0], A->m);
		printf(",\n");
		for (size_t i = 1; i < 3; i++) {
			printf("	");
			_print_ld(A->val[i], A->m);
			printf(",\n");
		}
		printf( "	...,\n");

		for (size_t i = A->n - 3; i < A->n - 1; i++) {
			printf("	");
			_print_ld(A->val[i], A->m);
			printf(",\n");
		}
		printf("	");
		_print_ld(A->val[A->n - 1], A->m);
	}
	printf("])\n");
}

/**
 * @brief Print imatrix to file (Python numpy library format). 
 * If row size or column size larger than 20, print brief representation
 * @param A integer matrix to print
 */
void fprint_imatrix(FILE *fp, imatrix *A) {

	#ifdef DEBUG
	fprintf(fp, "shape (%ld, %ld)\n", A->n, A->m);
	#endif

	fprintf(fp, "matrix([");
	if (A->n <= 20) {
		_fprint_ld(fp,A->val[0], A->m);
		fprintf(fp, ",\n");	
		for (size_t i = 1; i < A->n - 1; i++) {
			fprintf(fp, "	");
			_fprint_ld(fp,A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp, "	");
		_fprint_ld(fp,A->val[A->n - 1], A->m);
	} else {
		_fprint_ld(fp,A->val[0], A->m);
		fprintf(fp, ",\n");
		for (size_t i = 1; i < 3; i++) {
			fprintf(fp, "	");
			_fprint_ld(fp,A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp,  "	...,\n");

		for (size_t i = A->n - 3; i < A->n - 1; i++) {
			fprintf(fp, "	");
			_fprint_ld(fp,A->val[i], A->m);
			fprintf(fp, ",\n");
		}
		fprintf(fp, "	");
		_fprint_ld(fp,A->val[A->n - 1], A->m);
	}
	fprintf(fp, "])\n");
}

/**
 * @brief Read csv or dat file (one column file only, use for structure.dat)
 * @param filename name of csv or dat file to read
 * @return imatrix* 
 */
imatrix *read_icsv(char *filename) {

	FILE *fp;
	char line[1024];
	size_t n = 0, m = 1;
	size_t j = 0;
	int counter = 0;
    imatrix *output;

	fp = fopen(filename, "r");
 
    if (!fp) {
        printf("Can't open file\n");
        
        return NULL;
    } else {

		#ifdef DEBUG
		printf("\rRead file %s", filename);
		fflush(stdout);
		#endif
		
		while (fgets(line, sizeof(line), fp)) {
			n++;
		}

		/* Find number of separators to deduce number of columns */
		// TODO: Error handling. Inconsistent number of columns
		while(line[j] != '\0') {
			if(line[j] == ',')
				counter++;
			j++;
		}

		m = counter + 1;
				
		/* Reset pointer to beginning of file after reading number of
		 * lines */
		fseek(fp, 0, SEEK_SET);

		output = init_imatrix(n, m);
				
		for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++)
			    fscanf(fp, "%ld%*[,\n]", &output->val[i][j]);
		}
		
		fclose(fp);
		
		#ifdef DEBUG
		printf("\r[ok]Read file %s\n", filename);
		#endif

		return output;
	}
}

/**
 * @brief Print integer vector to file
 * @param v integer vector to save in file
 * @param filename name of file
 * @return int 0 if OK, 1 if cannot open file
 */
int ivtofile(ivector *v, const char *filename) {

	FILE *fp;
	size_t i = 0;

	fp = fopen(filename, "w");

	if (!fp) {
		printf("Cannot open file\n");

		return 1;
	} else {
		for (i = 0; i < v->n; i++)
			fprintf(fp, "%ld\n", v->val[i]);

		fclose(fp);

		return 0;
	}
}

/**
 * @brief Print integer matrix to file
 * @param v integer matrix to save in file
 * @param filename name of file
 * @return int 0 if OK, 1 if cannot open file
 */
int imtofile(imatrix *M, const char *filename) {

	FILE *fp;
	size_t i = 0, j = 0;

	fp = fopen(filename, "w");

	if (!fp) {
		printf("Cannot open file\n");

		return 1;
	} else {
		for (i = 0; i < M->n; i++) {
			for (j = 0; j < M->m - 1; j++)
				fprintf(fp, "%ld,", M->val[i][j]);
			fprintf(fp, "%ld\n", M->val[i][j]);
		}

		fclose(fp);

		return 0;
	}
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
 * @brief Compute matrix of absolute values
 * @param A integer matrix
 * @return imatrix* 
 */
imatrix *imabs(imatrix *A) {

	size_t i = 0, j = 0;
	imatrix *output;

	output = init_imatrix(A->n, A->m);

	for (i = 0; i < A->n; i++) {
		for (j = 0; j < A->m; j++) {
			output->val[i][j] = labs(A->val[i][j]);
		}
	}

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

	switch (p) {
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
 * @brief Compute p-distance between matrix A and matrix B
 * @param A integer matrix
 * @param B integer matrix
 * @param p 0 < p <= Inf such that we get ||A - B||_p (enter 0 for p = Inf)
 * @return double 
 */
double imdist(imatrix *A, imatrix *B, unsigned int p) {

	double dist = 0;
	imatrix *tmp, *C;

	/* Compute A - B first and then |A - B| */
	tmp = imsub(A, B);
	C = imabs(tmp);
	dist = imnorm(C, p);

	free_imatrix(tmp);
	free_imatrix(C);

	return dist;
}

/**
 * @brief Return integer vector of size n of increasing numbers, starting 0
 * @param n size of vecto
 * @return *ivector 
 */
ivector *iarange(size_t n) {

	ivector *v = init_ivector(n);

	for (size_t i = 0; i < n; i++)
        (v->val)[i] = i;

	return v;
}

/**
 * @brief return integer vector of ones of size 
 * @param n size of vector 
 * @return *ivector
 */
ivector *iones(size_t n) {

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
imatrix *ieye(size_t n) {

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
imatrix *idiag(ivector *v) {

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
long ivdot(ivector *v1, ivector *v2) {

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
 * @brief return vector substraction of two vectors
 * @param a integer vector a
 * @param b integer vector b
 * @return *ivector
 */
ivector *ivsub(ivector *a, ivector *b) {

    // Assert ivector have same size
	assert(a->n == b->n);

    ivector *c = init_ivector(a->n);

    for (size_t i = 0; i < c->n; i++)
        (c->val)[i] = (a->val)[i] - (b->val)[i];

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
 * @brief return matrix substraction of two integer matrices
 * @param A integer matrix A
 * @param B integer matrix B
 * @return *imatrix
 */
imatrix *imsub(imatrix *A, imatrix *B) {

    // Assert imatrix have same dimensions
	assert((A->n == B->n) && (A->m == B->m));

    imatrix *C = init_imatrix(A->n, A->m);

    for (size_t i = 0; i < C->n; i++) {
        for (size_t j = 0; j < C->m; j++) {
            (((C->val)[i]))[j] = (((A->val)[i]))[j] - (((B->val)[i]))[j];
        }
    }

    return C;
}

/**
 * @brief Compute coordinates of maximum element in integer matrix A
 * @param A integer matrix
 * @return size_t* 
 */
size_t *immax(imatrix *A) {

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
 * @brief Find element candidate in vector v and return vector of same size 
 * containing 1 if candidate is in position i and 0 otherwise
 * @param v integer vector to search
 * @param candidate element to search in vector
 * @return ivector* 
 */
ivector *ivfind(ivector *v, long candidate) {

	size_t k = 0, counter = 0;
	ivector *indices;

	indices = init_ivector(v->n);
	 
	for (k = 0; k < v->n; k++) {
		if (v->val[k] == candidate) {
			indices->val[k] = 1;
			counter++;		
		} else {
			indices->val[k] = 0;
		}
	}

	#ifdef DEBUG
	printf("Found %ld / %ld\n", counter, v->n);
	#endif
	
	return indices;
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
            (((C->val)[i]))[j] = c_ij;
            c_ij = 0;
        }
    }

    return C;
}

/**
 * @brief Sum coordinates of vector v between coordinate a and b 
 * (inclusive)
 * @param v integer vector to sum
 * @param a lower coordinate
 * @param b upper coordinate
 * @return long int
 */
long int isum(ivector *v, size_t a, size_t b) {

	// Assert coordinates are in range
	assert(((a >= 0) && (b < v->n)) && (a <= b));
	
	if (a == b)
		return 0;

	long int s = 0;
	
	for (long int i = a; i <= b; i++)
		s += v->val[i];
		
	return s;
}

/**
 * @brief Transpose integer matrix A
 * @param A integer matrix to transpose
 * @return imatrix*
 */
imatrix *itranspose(imatrix *A) {

	imatrix *B = init_imatrix(A->m, A->n);

	for (size_t i = 0; i < A->n; i++) {
		for (size_t j = 0; j < A->m; j++) {
			B->val[i][j] = A->val[j][i];
		}
	}

	return B;
}

/**
 * @brief Transform single-column integer matrix to an integer vector. Free 
 * matrix pointer.
 * @param A matrix to vectorize
 * @return *ivector
 */
ivector *ivectorize(imatrix *A) {

	ivector *output;

	// Check matrix has one column only
	assert(A->m == 1);

	output = init_ivector(A->n);

	for (size_t i = 0; i < A->n; i++)
		output->val[i] = A->val[i][0];

	free_imatrix(A);

	return output;
}