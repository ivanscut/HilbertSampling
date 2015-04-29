#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <sstream>
#include "generate_RQMC.h"
#include "HilbertSampling.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
//                                                                  //
//This code was used to carry out the results in He and Owen (2014).//
//                                                                  //
//////////////////////////////////////////////////////////////////////

extern "C" {

typedef void Fn(double *, int *, double *);//the prototype of integrands
void IndicatorFn(double *, int *, double *);
void sumFn(double *, int *, double *);
void cuspFn(double *, int *, double *);
string chooseIntegrand();
void oneSobolSimulation(int *, int *, int *, int *, double *, double *);
void batchSobolSimulation(int *, int *, int *, int *, double *, double *, double *);
void oneHilbertSimulation(int *,  int *, int *, double *, double *);
void batchHilbertSimulation(int *, int * , int *, double *, double *, double *);
void printResult(double *, double *, int );
void saveVector(double *, int , ofstream& );

/*========= the list of integrands==========*/
void IndicatorFn(double * cp, int *d, double * fx){
	double tmp = 0.0;
	double theta = (double) (*d)/2.0;
	for(int i=0; i< *d; i++)
		tmp += cp[i];
	fx[0] = tmp>theta? 1.0:0.0;
}
void sumFn(double * cp, int *d, double * fx){
	double tmp = 0.0;
	for(int i=0; i< *d; i++)
	 tmp += cp[i];
	fx[0] = tmp;
}
void cuspFn(double * cp, int *d, double * fx){
	double tot = 0.0;
	double tmp = (double) (*d)/2.0;
	for(int i=0; i< *d; i++)	
		tot += cp[i];	
	fx[0] = max(tot-tmp,0.0);
}
/*=========choose the integrand============*/
Fn  *myFn = NULL;
string chooseIntegrand()
{
	int opt;
	string name;
	cout << "Choose an integrand: smooth[1], cups[2], discontinuous[3]" << endl;
	cin >> opt;
	switch(opt)
	{
		case 1:
			myFn = sumFn;
			cout << "You are testing the smooth function..."<< endl;
			name = "smooth";
			break;
		case 2:
			myFn = cuspFn;
			cout << "You are testing the cusp function..."<< endl;
			name = "cusp";
			break;
		case 3:
			myFn = IndicatorFn;
			cout << "You are testing the discontinuous function..."<< endl;
			name = "discontinuous";
			break;
		default:
			myFn = IndicatorFn;
			cout << "You are testing the discontinuous function..."<< endl;
			name = "discontinuous";		
	}
	return name;
}

void oneSobolSimulation(int *sq, int *d, int * R, int * ns, double *ave, double *var){

	int N = (int) 1 << (*ns);
	double* p = new double[(*d)*N];
	int flag = 0;
	double tot;
	double m1 = 0.0, m2 = 0.0;
	Scrambled * os = Scrambled_Create(sq, d, &N) ;
	double fx;
	double * cp = new double [*d];
	for(int r = 0;r < (*R) ;r++){		
		Scrambled_Randomize(os);
		Scrambled_GetPoints(os, d, &N, p, &flag);
		tot = 0.0;
		for(int i = 0; i< N; ++i){			
			for(int j = 0; j< (*d); ++j){
				cp[j] = p[i*(*d)+j];
			}
			/*=========evaluate the integrand=========*/
			myFn(cp, d, &fx);
			tot += fx; 
		}
		tot /= N;
		m1 += tot;
		m2 += tot*tot;
		
	}
	*ave = (double) m1/(*R);
	*var = (double) (m2/(*R)-(*ave)*(*ave))*(*R)/((*R)-1);
	delete os;
}
void batchSobolSimulation(int *sq, int *d, int * R, int * ns, double *ave, double *var, double *cost){
	clock_t start_time,end_time;
	for(int i = 0; i <= (*ns); ++ i)
	{
		start_time = clock();	
		oneSobolSimulation(sq, d, R, &i,&ave[i], &var[i]);
		end_time = clock();	
		cost[i] = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;	
	}
}
void oneHilbertSimulation(int *d,  int *R, int *ns, double *ave, double *var){
	int N = (int) 1<< (*ns);
	double tot;
	int sq = 1, flag = 0;
	int order = (int) ceil((double) (*ns)/(*d));
	
	double m1 = 0.0, m2 = 0.0;
	int d0 = 1;
	double* vdc = new double[N];
	Scrambled * vdc_os = Scrambled_Create(&sq,&d0, &N) ;
	
	int order2  = (int) floor((double) (*ns)/(*d));
	HilbertGenerator os(order,order2,*d);	//the fast version
	//HilbertGenerator os(order,*d);
	
	double fx;
	double * cp = new double [*d];	
	for(int r = 0;r < (*R) ;r++){			
		Scrambled_Randomize(vdc_os);
		Scrambled_GetPoints(vdc_os, &d0, &N, vdc, &flag);
		tot = 0.0;
		for(int i = 0; i< N; ++i){					
			os.fastConvert(cp, vdc[i]);		
			//os.convert(cp, vdc[i]);					
			/*========evaluate the integrand=======*/
			myFn(cp, d, &fx);
			tot += fx;			
		}		
		tot /= N;
		m1 += tot;
		m2 += tot*tot;
	}
	*ave = (double) m1/(*R);
	*var = (double) (m2/(*R)-(*ave)*(*ave))*(*R)/((*R)-1);
	delete vdc_os;
}
void batchHilbertSimulation(int *d, int * R, int * ns, double *ave, double *var, double *cost){
	clock_t start_time,end_time;
	for(int i = 0; i <= (*ns); ++ i)
	{
		start_time = clock();	
		oneHilbertSimulation(d, R, &i, &ave[i], &var[i]);	
		end_time = clock();	
		cost[i] = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;
	}
}

void printResult(double *ave, double *var, int N){
	for(int i=0;i<=N;i++)
	{
		cout <<"log(n) = "<< i<< " Mean: " << ave[i]  << " Var: " << var[i] << endl;
	}
}
void saveVector(double * data, int N, ofstream& myfile){
    for(int i = 0;i < N;++i)
	{
		myfile << data[i] << " ";
	}	
    myfile << endl;
}
#define MSES //here choose the main to run

#ifdef RUNTIMES
int main(){

	int d, n, R, D;
	char filename[64];	
	sprintf(filename, "Runtime.txt");	
	ofstream myfile(filename);
	clock_t start_time,end_time;
	double t;
	int N;
	double * cp, *p, *vdc;
	int sq = 1;
	int flag = 0;
	cout << "Enter: log2(n), repetition, maximum dimension (1 to 31)" << endl;
	cin >> n >> R >> D;
	N = (int) 1<<n;
	for(d = 2; d <= D; ++d){			
		cp = new double [d];
		p  = new double [N*d];
		vdc = new double [N];
		myfile << d << " ";	
		cout << "You are testing d = " << d << "..."<< endl;
		cout << "Sobol' sampling starts..." << endl;		
		start_time = clock();
		Scrambled * os1 = Scrambled_Create(&sq, &d, &N) ;	
		for(int r = 0; r < R; ++r)
		{
			Scrambled_Randomize(os1);
			Scrambled_GetPoints(os1, &d, &N, p, &flag);
			for(int i = 0; i < N; ++i)
			{			
				for(int j = 0; j< d; ++j){
					cp[j] = p[i*d+j];
				}
			}
		}
		end_time = clock();		
		t = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;
		myfile << t << " ";
		cout << "t: " << t <<endl;
	
		cout << "Hilbert sampling starts..." << endl;		
		start_time = clock();
		int order = (int) ceil((double) n/d);
		int order2 = (int) floor((double) n/d);
		HilbertGenerator os(order,order2,d);//the fast version
		//HilbertGenerator os(order,d);
		int d0 = 1;
		Scrambled * vdc_os = Scrambled_Create(&sq,&d0, &N) ;
		for(int r = 0; r < R; ++r)
		{
			Scrambled_Randomize(vdc_os);
			Scrambled_GetPoints(vdc_os, &d0, &N, vdc, &flag);			
			for(int i = 0;i < N; ++i)
			{			
				os.fastConvert(cp,vdc[i]);
				//os.convert(cp,vdc[i]);
			}
		}
		end_time = clock();
		t = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;
		myfile << t << " ";
		cout << "t: " << t <<endl;		
		myfile << endl;
	}		
	return 0;
}
#endif

#ifdef MSES
int main(){	
	int d, N, R, sq = 1;	
	string fnname;
	stringstream ss("");	
	double t;
	while(1){
	fnname = chooseIntegrand();
	cout <<"Enter dimension (1 to 31), log2(n), repetition" <<endl;
	cin >> d >> N >> R;
	ss.str("");
	ss << fnname << "_d" << d << "_n" << N << "_R" << R << ".txt";
	
	ofstream myfile(ss.str().c_str());
	double *ave = new double[N+1];
	double *var = new double[N+1];
	double *count = new double[N+1];
	double *cost = new double[N+1];

	clock_t start_time,end_time;
	/*======test the Sobol' sampling========*/
/*	
	cout << "Sobol' sampling starts..." << endl;
	start_time = clock();
	batchSobolSimulation(&sq, &d, &R, &N, ave, var,cost);
	end_time = clock();	
	t = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;	
	cout<< "Run time is: "<< t <<"ms"<<endl;//total run time
	printResult(ave,var,N);	
	saveVector(var, N+1, myfile);
	saveVector(cost, N+1, myfile);
*/
	/*======test the Hilbert sampling========*/
	cout << "Hilbert sampling starts..." << endl;
	start_time = clock();
	batchHilbertSimulation(&d, &R, &N, ave, var,cost);
	end_time = clock();
	t = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000;	
	cout<< "Run time is: "<< t <<"ms"<<endl;//total run time
	printResult(ave,var,N);	
	saveVector(var, N+1, myfile);
	saveVector(cost, N+1, myfile);	
	
	/*======asking the next action======*/
	int continued;
	cout << "continue[1], stop[else]" << endl;
	cin>> continued;
	if(continued!=1)
		break;	
	}	
	system("pause");
	return 0;
}
#endif
}
