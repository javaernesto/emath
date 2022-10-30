#ifndef RANDOM_H
#define RANDOM_H

#include "fvector.h"
#include "ivector.h"


/* Function prototypes */

ivector *randsample(size_t n, int lower, int upper);

ivector *randsample_weight(size_t n, int lower, int upper, fvector *weights);

long randi_weight(int lower, int upper, fvector *weights);

void ishuffle(ivector *v);

#endif