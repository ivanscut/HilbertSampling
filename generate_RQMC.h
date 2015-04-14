#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include "DigitalNetsBase2.h"
#include "Drand48.h"

#define MINU 1e-10
typedef unsigned int  uint;

//////////////////////////////////////////////////////////////////////
//                                                                  //
// These are the interfaces with the QMC generators of            //
// Thomas Kollig and Alexander Keller. You can Download it through  //
// http://www.uni-kl.de/AG-Heinrich/SamplePack.html                 //
//                                                                  //
//////////////////////////////////////////////////////////////////////

bitvector rng32()
{
  return (bitvector)(drand48()*((double)BV_MAX+1.0));
}

extern "C" {
// interface for the constructor

DigitalNetGenerator *DigitalNetGenerator_Create(int *seq, int *d) 
{
	dgt gentype;

	if(*seq==1)
	{
		gentype=dgt_Sobol;
	}
	else if(*seq==2)
	{
		gentype=dgt_SpecialNiederreiter;
	}
	else if(*seq==3)
	{
		gentype=dgt_NiederreiterXing;
	}
    	return new DigitalNetGenerator(gentype, *d);
}

LinearScrambled * LinearScrambled_Create(int *seq, int *d) 
{
	dgt gentype;

	if(*seq==1)
	{
		gentype=dgt_Sobol;
	}
	else if(*seq==2)
	{
		gentype=dgt_SpecialNiederreiter;
	}
	else if(*seq==3)
	{
		gentype=dgt_NiederreiterXing;
	}
	else
	{
		gentype=dgt_ShiftNet;
	}
	
    	return new LinearScrambled(gentype, *d, rng32);
}
 

Scrambled * Scrambled_Create(int *seq, int *d, int *N) 
{
	dgt gentype;

	if(*seq==1)
	{
		gentype=dgt_Sobol;
	}
	else if(*seq==2)
	{
		gentype=dgt_SpecialNiederreiter;
	}
	else if(*seq==3)
	{
		gentype=dgt_NiederreiterXing;
	}
	else
	{
		gentype=dgt_ShiftNet;
	}
	
    return new Scrambled(gentype, *d, *N+1, rng32);
}


// interface to randomize points set

void LinearScrambled_Randomize(LinearScrambled *os) 
{
    os->Randomize();
}


void Scrambled_Randomize(Scrambled *os) 
{
    os->Randomize();
}

// when flag = 1, the coordinates are truncated in [MINU, 1-MINU]; otherwise, they are in [0,1].
void DigitalNetGenerator_GetPoints(DigitalNetGenerator *os, int* d, int* N, double *x, int *flag) 
{
	int i,j;
	
	for (i = 0; i < *N; i++)
	{
	  	for (j = 0; j < *d; j++)
		{
			if(*flag){
	    		x[i*(*d)+j] = max((*os)[j],MINU);
				x[i*(*d)+j] = min((*os)[j],1-MINU);
			}
			else
				x[i*(*d)+j] = (*os)[j];
		}
		++(*os);
	}

}

void LinearScrambled_GetPoints(LinearScrambled *os, int* d, int* N, double *x, int *flag)  
{
	int i,j;
	
	++(*os);
	for (i = 0; i < *N; i++)
	{
	  	for (j = 0; j < *d; j++)
		{
	    	if(*flag){
				x[i*(*d)+j] = max((*os)[j],MINU);
				x[i*(*d)+j] = min((*os)[j],1-MINU);
			}
			else
				x[i*(*d)+j] = (*os)[j];
		}
		++(*os);
	}

}


void Scrambled_GetPoints(Scrambled *os, int* d, int* N, double *x, int *flag) 
{
	int i,j;
	
	for (i = 0; i < *N; i++)
	{
	  	for (j = 0; j < *d; j++)
		{
	    	if(*flag){
	    		x[i*(*d)+j] = max((*os)[j],MINU);
				x[i*(*d)+j] = min((*os)[j],1-MINU);
			}
			else
				x[i*(*d)+j] = (*os)[j];
		}
		++(*os);
	}

}
}

