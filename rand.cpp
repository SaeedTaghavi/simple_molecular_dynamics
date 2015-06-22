#include "rand.hpp"

void randomize(unsigned int seed) {
	srandom(seed);
}

void randomize() {
	srandom(time(0) * getpid());
}

double drand(double inf, double sup) {
	return (sup - inf) * double(random()) / double(RAND_MAX) + inf;
}

int irand(int inf, int sup) {
	return int(floor(drand(inf, sup )));
}

