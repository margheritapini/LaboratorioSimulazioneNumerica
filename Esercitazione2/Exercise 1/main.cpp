/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department - Universita' degli Studi di Milano
  _/  _/_/    _/    _/       Exercises 01 
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
   
   
   // Esercizio 02.1.1: computing the integral sampling a uniform distribution in [0,1]
   
   ofstream output1;
   output1.open("integral_values_unif.dat");
   
   int M = 10000000; // Total number of throws
   int N = 100; // Number of blocks
   int L = int(M/N); // Number of throws in each block
   
   vector<double> integral_value(N,0);
   vector<double> integral_value2(N,0);
   double sum_integral = 0;
   double x = 0; // Variabile della funzione 
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	sum_integral = 0;
      	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		x = rnd.Rannyu();
   		sum_integral += M_PI/2.*cos(M_PI/2*x);
     	}
   	
   	integral_value[i] = sum_integral/L;
   	integral_value2[i] = integral_value[i]*integral_value[i];   	
      }
      
   double sum_prog = 0;
   double su2_prog = 0;
   double err_prog = 0;
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog = 0;
   	su2_prog = 0;
   	err_prog = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog += integral_value[j];
   		su2_prog += integral_value2[j];
   	}
   	
   	sum_prog /= (i+1); // Cumulative average
    	su2_prog /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog = 0;
    		}
    	else {
    		err_prog = sqrt((su2_prog - sum_prog*sum_prog)/i); 
    		}
    	
    	output1 << L*i << " " << sum_prog << " " << err_prog << " " << endl;
	
	}
	
	output1.close();
	
	// Esercizio 02.1.2: computing the integral using importance sampling in [0,1]
	
	ofstream output2;
   output2.open("integral_values_importance.dat");
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	sum_integral = 0;
      	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		x = rnd.my_prob();
   		sum_integral += M_PI/2.*cos(M_PI/2*x)/(2-2*x);
     	}
   	
   	integral_value[i] = sum_integral/L;
   	integral_value2[i] = integral_value[i]*integral_value[i];   	
      }
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog = 0;
   	su2_prog = 0;
   	err_prog = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog += integral_value[j];
   		su2_prog += integral_value2[j];
   	}
   	
   	sum_prog /= (i+1); // Cumulative average
    	su2_prog /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog = 0;
    		}
    	else {
    		err_prog = sqrt((su2_prog - sum_prog*sum_prog)/i); 
    		}
    	
    	output2 << L*i << " " << sum_prog << " " << err_prog << " " << endl;
	
	}
	
	output2.close();
   
   
	return 0;
}

