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

double funct_chi2 (double, double, double);

double funct_chi2 (double n_i, double n_tot, double M) {
	return (n_i-n_tot/M)*(n_i-n_tot/M)/(n_tot/M);
}

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
   
   // 1 - 2
   
   ofstream output1;
   ofstream output2;
   
   output1.open("N_integral_values.dat");
   output2.open("N_sigma2_values.dat");
   
   int M = 10000000; // Total number of throws
   int N = 100; // Number of blocks
   int L = int(M/N); // Number of throws in each block
   
   
   // For each i-th block:
   // - I compute the integral value r_i (obtained as the mean value of the L extracted random numbers)
   // - I compute the double square of that value, which I will use to compute the error for each iteration
   
   double ave[N] = {}; // Array for the integral value r_i
   double av2[N] = {}; // Array for r_i*r_i
   
   double ave_sigma[N] = {}; // Array for the integral value r_i
   double av2_sigma[N] = {}; // Array for r_i*r_i
      
   double sum = 0;
   double sum_sigma = 0;
   
   double rand_numb = 0;
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   
   	sum = 0; // Cumulative variable to contain the sums of the extracted random numbers
   	sum_sigma=0;
   	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		
   		rand_numb = rnd.Rannyu();
   		sum += rand_numb;
   		sum_sigma += (rand_numb - 0.5)*(rand_numb - 0.5);
   		}
   	
   	ave[i] = sum/L; // r_i: integral value for the i-th block
   	av2[i] = ave[i]*ave[i]; // r_i*r_i
   	
   	ave_sigma[i] = sum_sigma/L; // r_i: integral value for the i-th block
   	av2_sigma[i] = ave_sigma[i]*ave_sigma[i]; // r_i*r_i
   	   
   } 
   
   
   // Now I have to compute the value of the integral as the number of iterations increases. For L throws, the value of the integral is r_0. For 2L throws, the value of the integral is (r_0+r_1)/2. For N*L throws, the value of the integral is (r_0+r_1+...+r_N-1)/N.
   
   
   double sum_prog = 0;
   double su2_prog = 0;
   double err_prog = 0;
   
   double sum_prog_sigma = 0;
   double su2_prog_sigma = 0;
   double err_prog_sigma = 0;
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   
   	sum_prog = 0;
   	su2_prog = 0;
   	err_prog = 0;
   	
   	sum_prog_sigma = 0;
   	su2_prog_sigma = 0;
   	err_prog_sigma = 0;
   	
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog += ave[j];
   		su2_prog += av2[j];
   		
   		sum_prog_sigma += ave_sigma[j];
   		su2_prog_sigma += av2_sigma[j];
   	}
   	
   	sum_prog /= (i+1); // Cumulative average
    	su2_prog /= (i+1); // Cumulative square average
    	
    	sum_prog_sigma /= (i+1); // Cumulative average
    	su2_prog_sigma /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog = 0;
    		err_prog_sigma = 0;}
    	else {
    		err_prog = sqrt((su2_prog - sum_prog*sum_prog)/i); 
    		err_prog_sigma = sqrt((su2_prog_sigma - sum_prog_sigma*sum_prog_sigma)/i); 
    		}
    		
    	output1 << i << " " << sum_prog << " " << err_prog << " " << endl;
    	output2 << i << " " << sum_prog_sigma << " " << err_prog_sigma << " " << endl;
   
   }
   
  
  output1.close(); 
  output2.close();
  
  // 3
  
  ofstream output3;
  output3.open("chi2_values.dat");
  
  double subintervals = 100; // Number of sub-intervals of [0,1]
  double n = 10000; // Number of throws for each iteration
  double iteration = 1000;
  
  vector<int> counts(subintervals, 0); // Vector to store the counts of numbers in each subinterval
  double chi2 = 0;
  
  int subinterval_index = 0;
  
  for (int j=0; j<iteration; j++) { // I repeat the test 100 times
  	
  	chi2 = 0; // Re-inizialization of the variables at every iteration
  	fill(counts.begin(), counts.end(),0);
  	
  	for (int i=0; i<n; i++) {
  		subinterval_index = static_cast<int>(rnd.Rannyu() * subintervals); 
  		counts[subinterval_index]++;
  	}
  	
  	// chi2 computation
  	
  	for (int i = 0; i<subintervals; i++) {
  		chi2 += funct_chi2 (counts[i], n, subintervals);
  	}
  	
  	output3 << j << " " << chi2 << endl;
  }
  
  output3.close();
  
  return 0;
}

