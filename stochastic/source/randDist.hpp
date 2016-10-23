

#ifndef RAND_DIST_H
#define RAND_DIST_H

#include <math.h>

// uniform distribution (returns a double from 0.0-1.0)
inline double unif_dist () {
	return (rand()+1) / (((double)RAND_MAX) + 2);
}

// binomial distribution (returns an integer from 0-n)
inline int bino_dist (int n, double p) {
	int x = 0;
	for (int i = 0; i < n; i++) {
		if (unif_dist() < p) {
			x++;
		}
	}
	return x;
}

// exponential distribution (returns a double based on mean)
inline double expo_dist (double mean) {
	return -log(1.0 - unif_dist()) / mean;
}

// Poisson distribution (returns an integer based on mean)
inline int pois_dist (double mean) {
	double L = pow(M_E, -mean);
	int k = 0;
	double p = 1;
	do {
		k++;
		p *= unif_dist();
	} while (p > L);
	return k - 1;
}

// Pk distribution, which is just a logarithmic uniform distribution
inline double pk_dist () {
	return log(1.0 / unif_dist());
}

#endif

