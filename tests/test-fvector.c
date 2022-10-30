#include <stdio.h>

#include "../include/fvector.h"

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