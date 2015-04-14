#ifndef _HILBERSAMPLING_H_
#define _HILBERSAMPLING_H_

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Version 1.0 (2015-04-112)                                               //
// Copyright (C) 2015 by Zhijian He(hezhijian87@gmail.com)                 //
//                       Art B. Owen(owen@stanford.edu)                    //
//                                                                         //
// All rights reserved. You may not distribute this software, in whole or  //
// in part, especially not as part of any commercial product, without the  //
// express consent of the authors.                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// This program uses the Hibert Mapping Algorithms in "Calculation         //
// of Mappings Between One and n-dimensional Values Using the Hilbert      //
// Space-filling Curve" by J. K. Lawder.                                   //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "Drand48.h"
typedef unsigned long      bitvector;
typedef unsigned int       uint;
typedef unsigned long long ulonglong;

typedef struct
{
    double*    value;
    ulonglong  xJ;
	ulonglong  W;
	bool       tag;
}HilbertData;

class HilbertGenerator
{
	private:				
		ulonglong   *rho;//derived-key 
		int         ORDER; //order of the Hilbert curve
		int         DIM; //dimension		
		int         resORDER; //order of the Hilbert curve to be stored
		uint        calc_J(ulonglong);
		ulonglong   calc_T(ulonglong);
		ulonglong   calc_tS_tT(ulonglong,ulonglong);
		void        double2DerivedKey(double); //convert decimal to derived-key 
		HilbertData *data; //which is used to store the Hilbert data 
		uint        ndata; // number of data
	public:
		HilbertGenerator(const int, const int, const int);		
		HilbertGenerator(const int, const int);		
		void convert(double *, double);
		void fastConvert(double *, double);
		~HilbertGenerator();
};

/*===========================================================*/
/* fast algorithm constructor    */
/*===========================================================*/
HilbertGenerator::HilbertGenerator(const int order, const int order2, const int dim)
{	
	ORDER = max(order,1);
	DIM = dim;
	resORDER = order2;
	rho = new ulonglong[order];
	if(resORDER > 0)
	{
		ndata = 1 << (resORDER*dim);
		data = new HilbertData[ndata];
		//initialization
		for(int i = 0; i < ndata; ++i)
		{
			data[i].value = new double[dim];
			data[i].xJ = 0;
			data[i].W = 0;
			data[i].tag = false;//unused marker
		}
	}else
		data = NULL;
}


/*===========================================================*/
/* the usual constructor  */
/*===========================================================*/
HilbertGenerator::HilbertGenerator(const int order, const int dim)
{	
	ORDER = max(order,1);
	DIM = dim;
	rho = new ulonglong[order];	
	data = NULL;
}
/*===========================================================*/
/* destructor */
/*===========================================================*/
HilbertGenerator::~HilbertGenerator()
{ 
  delete []rho;  
  if(data)
	delete []data;
}

/*===========================================================*/
/* convert double to derived-key */
/*===========================================================*/
void HilbertGenerator:: double2DerivedKey(double input)
{
	uint tmp = 1 << DIM;
	for(uint i = 0; i < ORDER; ++i)
	{
		input *= tmp;
		rho [i] = (ulonglong) input;
		input -= rho [i];
	}
}

/*===========================================================*/
/* calc_J */
/*===========================================================*/
uint HilbertGenerator::calc_J (ulonglong P)
{
	uint i;
	uint J;
	J = DIM;
	for (i = 1; i < DIM; i++)
		if ((P >> i & 1) == (P & 1))
			continue;
		else
			break;
	if (i != DIM)
	J -= i;
	return J;
}
/*===========================================================*/
/* calc_T */
/*===========================================================*/
ulonglong HilbertGenerator::calc_T (ulonglong P)
{
	if (P < 3)
		return 0;
	if (P % 2)
		return (P - 1) ^ (P - 1) / 2;
	return (P - 2) ^ (P - 2) / 2;
}
/*===========================================================*/
/* calc_tS_tT */
/*===========================================================*/
ulonglong HilbertGenerator::calc_tS_tT(ulonglong xJ, ulonglong val)
{
	ulonglong retval, temp1, temp2;
	retval = val;
	if (xJ % DIM != 0)
	{
		temp1 = val >> xJ % DIM;
		temp2 = val << DIM - xJ % DIM;
		retval = temp1 | temp2;
		retval &= ((ulonglong)1 << DIM) - 1;
	}
	return retval;
}
/*=========================================================================*/
/* Key function: map one dimension to DIM dimensions through Hilbert curve */
/*=========================================================================*/
void HilbertGenerator::convert(double * pt, double input)
{
	//get the derived-key through the input	
	double2DerivedKey(input);		
	
	ulonglong A, W, S, tS, T, tT, J, P, xJ;
	double mask;
	int i, j, b;
	//initial step			
	xJ = 0;
	W = 0;		
	for(b = 0; b < DIM; ++b)
		pt[b] = 0;
	for(i = 0, mask = 0.5; i < ORDER; ++i)
	{
			P = rho[i];
			S = P^P/2;
			tS = calc_tS_tT(xJ, S);
			A = W^tS;
			/*--- distrib bits to coords ---*/
			for (j = DIM - 1; A > 0,j>=0; A >>=1, j--)
				if (A & 1)
					pt[j] += mask;
			T = calc_T(P);
			tT = calc_tS_tT(xJ, T);
			W ^= tT;
			J = calc_J(P);
			xJ += J - 1;
			mask /= 2;
	}
	//add the IID random offsets
	for(b = 0; b < DIM; ++b)
	{
		pt[b] += drand48()/(1<<ORDER);
	}
}


/*=========================================================================*/
/* Key function (fast version): map one dimension to DIM dimensions        */
/*=========================================================================*/
void HilbertGenerator::fastConvert(double * pt, double input)
{
	if(resORDER <= 0)//that redirects to the usual case
	{
		convert(pt, input);
		return;
	}
	//get the derived-key through the input
	double2DerivedKey(input);
	
	ulonglong A, W, S, tS, T, tT, J, P, xJ;
	double mask;
	int i, j, b;
	uint id = (uint) (ndata*input);
		
	if(data[id].tag)//if the data was recorded.
	{
		//recover the first resORDER bits
		mask = pow(2.0,-resORDER);
		xJ = data[id].xJ;
		W  = data[id].W;
		for(b = 0; b < DIM; ++b)
			pt[b] = data[id].value[b];
		
		//produce the remaining bits.
		for (i = resORDER, mask /= 2; i < ORDER; ++i)
		{
			P = rho[i];
			S = P^P/2;		
			tS = calc_tS_tT(xJ, S);			
			A = W^tS;
			/*--- distrib bits to coords ---*/
			for (j = DIM - 1; A > 0,j>=0; A >>=1, j--)
				if (A & 1)
					pt[j] += mask;			
			T = calc_T(P);
			tT = calc_tS_tT(xJ, T);
			W ^= tT;
			J = calc_J(P);
			xJ += J - 1;
			mask /= 2;			
		}
	}else
	{				
		xJ = 0;
		W = 0;		
		for(b = 0; b < DIM; ++b)
			pt[b] = 0;
		for(i = 0, mask = 0.5; i < ORDER; ++i)
		{
			P = rho[i];
			S = P^P/2;
			tS = calc_tS_tT(xJ, S);
			A = W^tS;
			/*--- distrib bits to coords ---*/
			for (j = DIM - 1; A > 0,j>=0; A >>=1, j--)
				if (A & 1)
					pt[j] += mask;
			T = calc_T(P);
			tT = calc_tS_tT(xJ, T);
			W ^= tT;
			J = calc_J(P);
			xJ += J - 1;
			mask /= 2;
			if(resORDER == i+1)//store the first resORDER bits
			{
				data[id].tag = true;
				for(b = 0; b < DIM; ++b)
					data[id].value[b] = pt[b];
				data[id].xJ  = xJ;
				data[id].W = W;
			}
		}
	
	}
	//add the IID random offsets
	for(b = 0; b < DIM; ++b)
	{
		pt[b] += drand48()/(1<<ORDER);
	}
}

#endif
