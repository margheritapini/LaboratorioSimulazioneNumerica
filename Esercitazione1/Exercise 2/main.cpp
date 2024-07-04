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


double compute_mean_value_linear (int N, double S_N, Random& random) {
	S_N = 0;
	
	for (int i=0; i<N; i++) {
		S_N += random.Rannyu();}
		
	return S_N/N;
}

double compute_mean_value_exp (int N, double S_N, Random& random, double lambda) {
	S_N = 0;
	
	for (int i=0; i<N; i++) {
		S_N += random.Exp(lambda);}
		
	return S_N/N;
}

double compute_mean_value_lorentz (int N, double S_N, Random& random, double mu, double Gamma) {
	S_N = 0;
	
	for (int i=0; i<N; i++) {
		S_N += random.CauchyLorentz(mu, Gamma);}
		
	return S_N/N;
}


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
   
   int iteration = 10000;
   double S_N = 0; // Mean value
   double lambda = 1; // Exponential parameter
   double mu = 0; // Lorentzian parameter
   double Gamma = 1; // Lorentzian parameter
  
	ofstream output1;
	ofstream output2;
	ofstream output3;
	ofstream output4;
	
	output1.open("N1.dat");
	output2.open("N2.dat");
	output3.open("N10.dat");
	output4.open("N100.dat");
	
  
	for (int i=0; i<iteration; i++) {
  		output1 << compute_mean_value_linear (1, S_N, rnd) << " " << compute_mean_value_exp (1, S_N, rnd, lambda) << " " << compute_mean_value_lorentz (1, S_N, rnd, mu, Gamma) << endl;
  		output2 << compute_mean_value_linear (2, S_N, rnd) << " " << compute_mean_value_exp (2, S_N, rnd, lambda) << " " << compute_mean_value_lorentz (2, S_N, rnd, mu, Gamma) << endl;
  		output3 << compute_mean_value_linear (10, S_N, rnd) << " " << compute_mean_value_exp (10, S_N, rnd, lambda) << " " << compute_mean_value_lorentz (10, S_N, rnd, mu, Gamma) << endl;
  		output4 << compute_mean_value_linear (100, S_N, rnd) << " " << compute_mean_value_exp (100, S_N, rnd, lambda) << " " << compute_mean_value_lorentz (100, S_N, rnd, mu, Gamma) << endl;  		
	}
	
	output1.close();
	output2.close();
	output3.close();
	output4.close();
   
	return 0;
}

