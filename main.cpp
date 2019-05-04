#ifndef __MAIN_CPP
#define __MAIN_CPP

#include "black_scholes.h"
#include "interval_bisection.h"
#include "newton_raphson.h"
#include "brents_method.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main() {

  	srand(time(NULL));

  	double S = 100.0;
  	double K = 100.0;
  	double C_M = 10.5;
  	double r = 0.06;   
	double T = 1.0;  

  	double low_vol = 0.1;
	double high_vol = 0.30;
	double episilon = 0.001;
	double init = 0.3;
  	
  	const clock_t begin_time_ib = clock();

  	for(int i=0;i<100000;i++){
		S += pow(-1,i)*10;  
	  	K = S;  
	  	C_M += 0.1*pow(-1,i);
	  
	  	//Black-Scholes functor
	  	//cout<<"Before BSC "<<i<<endl;
	  	BlackScholesCall bsc(S, K, r, T);
	  	//cout<<"After BSC "<<i<<endl;
	  	// Calculate the implied volatility
	  	//cout<<"Before IBS "<<i<<endl;
	  	
	  	double sigma_ib = interval_bisection(C_M, low_vol, high_vol, episilon, bsc);
	  	
	  	//cout<<"After IBS "<<i<<endl;
	  	//cout<<"\nSpot and Strike Price: "<<S<<", "<<K<<endl;
	  	//cout<<"Option premium: "<<C_M<<endl;
	  	//cout << "Implied Vol using Interval Bisection: " << sigma_ib << std::endl;  	
  	}
  	
  	cout<<"Time Taken using Interval Bisection: "<<float(clock()-begin_time_ib)/(CLOCKS_PER_SEC*100000)<<endl;

  	const clock_t begin_time_ns = clock();

  	for(int i=0;i<100000;i++){
		S += pow(-1,i)*10;  
	  	K = S;  
	  	C_M += 0.1*pow(-1,i);
	  
	  	//Black-Scholes functor
	  	//cout<<"Before BSC "<<i<<endl;
	  	BlackScholesCall bsc(S, K, r, T);
  		double sigma_ns = newton_raphson<BlackScholesCall, 
	  								  &BlackScholesCall::option_price,
	  								  &BlackScholesCall::option_vega>(C_M, init, episilon, bsc); 
	  	
		init = sigma_ns;
		//cout << "Implied Vol using Newton-Raphson: " << sigma_ns << std::endl;
	}	
	cout<<"Time Taken using Newton-Raphson: "<<float(clock()-begin_time_ns)/(CLOCKS_PER_SEC*100000)<<endl;

	double sigma_bm = 0.201318;

	const clock_t begin_time_bm = clock();

  	for(int i=0;i<100000;i++){
		S += pow(-1,i)*10;  
	  	K = S;  
	  	C_M += 0.1*pow(-1,i);
	  
	  	//Black-Scholes functor
	  	//cout<<"Before BSC "<<i<<endl;
	  	//BlackScholesCall bsc(S, K, r, T);
  		sigma_bm = brents_method(S, K, r, T, C_M, low_vol, high_vol, episilon);

  		//low_vol = 0.90*sigma_bm;
  		//high_vol = 1.10*sigma_bm;
	  	
		//cout << "Implied Vol using Newton-Raphson: " << sigma_ns << std::endl;
	}	
	cout<<"Time Taken using Brent's Method: "<<float(clock()-begin_time_bm)/(CLOCKS_PER_SEC*100000)<<endl;


  	return 0;
}

#endif