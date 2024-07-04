/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department - Universita' degli Studi di Milano
  _/  _/_/    _/    _/       Exercises 03 
 _/    _/       _/ _/       
_/    _/  _/_/_/  _/_/_/_/ 
*****************************************************************
*****************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"


using namespace std;


int main (int argc, char *argv[]){

   Random rnd;
   
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();
   
   // Esercizio 03.1: Computing at time t=0 the European call-option price and put-option price by sampling directly the final asset price S(T) for a GMB(r, sigma2)
   
   ofstream output1;
   output1.open("call_put_direct.dat");
   
   int M = 1000000; // Total number of throws
   int N = 100; // Number of blocks
   int L = int(M/N); // Number of throws in each block
   
   double sum_call = 0;
   double sum_put = 0;
   
   double Z = 0; // Auxiliar randomic variable
   
   double S0 = 100; // Asset price at t=0
   double T = 1; // Delivery time
   double K = 100; // Strike price
   double r = 0.1; // Risk-free interest rate
   double sigma = 0.25; // Volatility
   double S = 0; // Price
   
   vector<double> call_value (N,0);
   vector<double> call_value2 (N,0);
   vector<double> put_value (N,0);
   vector<double> put_value2 (N,0);
   
	for (int i=0; i<N; i++) { // Iteration over the blocks
   	sum_call = 0;
   	sum_put = 0;
      	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		Z = rnd.Gauss(0,1);
   		
   		S = S0*exp((r-sigma*sigma*0.5)*T+sigma*Z*sqrt(T));
   		
   		sum_call += exp(-r*T)*max(0., S-K);
   		sum_put += exp(-r*T)*max(0., K-S);
     	}
     	
        	
   	call_value[i] = sum_call/L;
   	put_value[i] = sum_put/L;
   	call_value2[i] = call_value[i]*call_value[i];
   	put_value2[i] = put_value[i]*put_value[i];      	
      }
      
   double sum_prog_call = 0;
   double sum_prog_put = 0;
   double su2_prog_call = 0;
   double su2_prog_put = 0;
   double err_prog_call = 0;
   double err_prog_put = 0;
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog_call = 0;
   	su2_prog_call = 0;
   	err_prog_call = 0;
   	
   	sum_prog_put = 0;
   	su2_prog_put = 0;
   	err_prog_put = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog_call += call_value[j];
   		su2_prog_call += call_value2[j];
   		
   		sum_prog_put += put_value[j];
   		su2_prog_put += put_value2[j];
   	}
   	
   	sum_prog_call /= (i+1); // Cumulative average
    	su2_prog_call /= (i+1); // Cumulative square average
    	
    	sum_prog_put /= (i+1); // Cumulative average
    	su2_prog_put /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog_call = 0;
    		err_prog_put = 0;
    		}
    	else {
    		err_prog_call = sqrt((su2_prog_call - sum_prog_call*sum_prog_call)/i);
    		err_prog_put = sqrt((su2_prog_put - sum_prog_put*sum_prog_put)/i); 
    		}
    	
    	output1 << L*i << " " << sum_prog_call << " " << err_prog_call << " " << sum_prog_put << " " << err_prog_put << " " << endl;
	
	}
	
	output1.close();
	
	// Esercizio 03.1: Computing at time t=0 the European call-option price and put-option price by sampling the discretized GMB(r, sigma2) path of the asset price
   
   ofstream output2;
   output2.open("call_put_discretized.dat");
   
   double time_intervals = 100;
   double ti = T/time_intervals;
   
    
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	sum_call = 0;
   	sum_put = 0;
      	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		S = 100;
   		
   		for (int k=0; k<time_intervals; k++) {
   			Z = rnd.Gauss(0,1);
   			S *= exp((r-sigma*sigma*0.5)*ti+sigma*Z*sqrt(ti));
   			}
   	   		
   		sum_call += exp(-r*T)*std::max(0., S-K);
   		sum_put += exp(-r*T)*std::max(0., K-S);
   		
     	}
   	
   	call_value[i] = sum_call/L;
   	put_value[i] = sum_put/L;
   	call_value2[i] = call_value[i]*call_value[i];
   	put_value2[i] = put_value[i]*put_value[i];      	
      }
      
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog_call = 0;
   	su2_prog_call = 0;
   	err_prog_call = 0;
   	
   	sum_prog_put = 0;
   	su2_prog_put = 0;
   	err_prog_put = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog_call += call_value[j];
   		su2_prog_call += call_value2[j];
   		
   		sum_prog_put += put_value[j];
   		su2_prog_put += put_value2[j];
   	}
   	
   	sum_prog_call /= (i+1); // Cumulative average
    	su2_prog_call /= (i+1); // Cumulative square average
    	
    	sum_prog_put /= (i+1); // Cumulative average
    	su2_prog_put /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog_call = 0;
    		err_prog_put = 0;
    		}
    	else {
    		err_prog_call = sqrt((su2_prog_call - sum_prog_call*sum_prog_call)/i);
    		err_prog_put = sqrt((su2_prog_put - sum_prog_put*sum_prog_put)/i); 
    		}
    	
    	output2 << L*i << " " << sum_prog_call << " " << err_prog_call << " " << sum_prog_put << " " << err_prog_put << " " << endl;
	
	}
   
   output2.close();
   
	return 0;
}

