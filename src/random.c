#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "../include/random.h"


/**
 * @brief Generate random integer vector of size n with uniform distribution, 
 * with values between lower and upper
 * @param n size of random vector
 * @param lower lower bound for values
 * @param upper upper bound for values
 * @return ivector* 
 */
ivector *randsample(size_t n, int lower, int upper) {
	
	double r = 0.0;
	int l = 0, h = 0, m = 0;
	
	ivector *output = init_ivector(n);
	
	for (size_t i = 0; i < n; i++) {
		// Random number between 0 and 1
        r = rand() / ((double) RAND_MAX);
        l = lower;
        h = upper;
        
        while (h > l) {
			m = (h + l) / 2;
			if (r > 1.0 / n) 		/* Probability for uniform distribution */
				l = m + 1;
			else
				h = m;
		}
		
		output->val[i] = l;
    }
	
	return output;
}

/**
 * @brief Generate random integer vector of size n with weights, with values
 * between lower and upper
 * @param n size of random vector
 * @param lower lower bound for values
 * @param upper upper bound for values
 * @param weights double vector containing weights
 * @return ivector* 
 */
ivector *randsample_weight(size_t n, int lower, int upper, fvector *weights) {
	
	double r = 0.0;
	int l = 0, h = 0, m = 0;
	
	// Check if [lower, upper] has same size as weights
	assert((upper - lower + 1) == weights->n);
	
	// Compute cummulative density function (assume sum = 1)
	fvector *cdf = init_fvector(weights->n);
	
	cdf->val[0] = weights->val[0];
	for (int i = 1; i < cdf->n; i++) {
		cdf->val[i] = cdf->val[i - 1] + weights->val[i];
	}
	
	ivector *output = init_ivector(n);
	
	for (size_t i = 0; i < n; i++) {
		// Random number between 0 and 1
        r = rand() / ((double) RAND_MAX);
        l = lower;
        h = upper;
        
        while (h > l) {
			m = (h + l) / 2;
			if (r > cdf->val[m-1])
				l = m + 1;
			else
				h = m;
		}
		
		output->val[i] = l;
    }

	free_fvector(cdf);
	
	return output;
}

long randi_weight(int lower, int upper, fvector *weights) {

	double r = 0.0;
	int l = 0, h = 0, m = 0;
	
	// Check if [lower, upper] has same size as weights
	assert((upper - lower + 1) == weights->n);
	
	// Compute cummulative density function (assume sum = 1)
	fvector *cdf = init_fvector(weights->n);
	
	cdf->val[0] = weights->val[0];
	for (int i = 1; i < cdf->n; i++)
		cdf->val[i] = cdf->val[i - 1] + weights->val[i];

	// Random number between 0 and 1
	r = rand() / ((double) RAND_MAX);
	l = lower;
	h = upper;
	
	while (h > l) {
		m = (h + l) / 2;
		if (r > cdf->val[m-1])
			l = m + 1;
		else
			h = m;
	}
		
	free_fvector(cdf);
	
	return l;
}

/**
 * @brief Random shuffle of integer vector v
 * @param v vector to shuffle
 */
void ishuffle(ivector *v) {

    if (v->n > 1) {
        size_t i;
        for (i = 0; i < v->n - 1; i++) {
          size_t j = i + rand() / (RAND_MAX / (v->n - i) + 1);
          int t = v->val[j];
          v->val[j] = v->val[i];
          v->val[i] = t;
        }
    }
}


ivector *randcycle(size_t n) {

	/* Assert cycle larger than id_1 cycle */
	assert(n > 1);

	size_t i = 0;
	ivector *cycle;

	cycle = init_ivector(n);

	for (i = 0; i < n - 1; i++) {
		size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
		cycle->val[i] = j;
	}

	return cycle;
}