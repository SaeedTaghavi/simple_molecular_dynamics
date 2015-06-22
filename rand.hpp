#ifndef _AZAR_UTILS_H_
#define _AZAR_UTILS_H_

#include <iostream>
#include <cstdlib>  
#include <cmath>
#include <ctime>
#include <unistd.h>

/***
 * INIT SEED
 */
void randomize();
void randomize(unsigned int) ;

/***
 * Return a random number in [ min , max ]
 */
double drand(double min = 0.0e0, double max = 1.0e0);
int irand(int min, int max) ;

#endif
