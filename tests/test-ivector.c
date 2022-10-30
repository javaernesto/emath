#include <stdio.h>

#include "../include/ivector.h"

int main() {

    // TESTING
	long v1[5] = {8, 2, 3, 4, 5};
	long v2[5] = {3, 0, 4, 0, 0};

	ivector *a = to_ivector(v1, 5);
	ivector *b = to_ivector(v2, 5);

    ivector *c = ivadd(a, b);

	printf("vdot %ld \n", ivdot(a, b));

    for (size_t i = 0; i < c->n; i++) {
        printf("%ld ", (c->val)[i]);
    }
    printf("\n");

	imatrix *mat = idiag(a);
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