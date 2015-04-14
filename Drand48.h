#ifndef DRAND48_H
#define DRAND48_H

#include <stdlib.h>
#include <functional>
#define m 0x100000000LL
#define c 0xB16
#define a 0x5DEECE66DLL
#define LB 1e-10
static unsigned long long seed = 1000;

double drand48(void)
{
	seed = (a * seed + c) & 0xFFFFFFFFFFFFLL;
	unsigned int x = seed >> 16;
    return 	((double)x / (double)m);
	
}

double getRand(int flag)
{
	double rng = drand48();
	if(flag)
	{
		rng = (rng==0)? LB: rng;
		rng = (rng==1)? 1-LB: rng;
	}
	return rng;
}

void srand48(unsigned int i)
{
    seed  = (((long long int)i) << 16) | rand();
}

#endif