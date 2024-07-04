/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department - Universita' degli Studi di Milano
  _/  _/_/    _/    _/       Exercises 02
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
   
   // Esercizio 02.2.1: random walk on a cubic lattice
   
   ofstream output1;
   output1.open("random_walk_lattice_blocks.dat");
   
   ofstream output2;
   output2.open("random_walk_lattice.dat");
   
   int M = 10000; // Total number of throws
   int N = 100; // Number of blocks
   int L = int(M/N); // Number of throws in each block
   
   int N_steps_rw = 100; // Number of steps random walk
   
   vector<double> position(3,0);
   double y = 0; // Auxiliar random number
   double a = 1; // step
   
   vector<double> random_walk_sum(N,0);
   vector<vector<double>> random_walk(N,vector<double>(N,0));
   vector<vector<double>> random_walk2(N,vector<double>(N,0));
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	fill(random_walk_sum.begin(), random_walk_sum.end(), 0);
   	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		fill(position.begin(), position.end(), 0);
   		
   		for (int k=0; k<N_steps_rw; k++) {
   			y = rnd.Rannyu(0,6);
   			
   			if (y<=1) position[0]+= a;
   			if (y>1 and y<=2) position[0] -= a;
   			if (y>2 and y<=3) position[1] += a;
   			if (y>3 and y<=4) position[1] -= a;
   			if (y>4 and y<=5) position[2] += a;
   			if (y>5) position[2] -= a;
   			
   			random_walk_sum[k] += position[0]*position[0]+position[1]*position[1]+position[2]*position[2]; 
   		}
   	}
   	
   	for (int k=0; k<100; k++) {
   	
   		random_walk[i][k] = sqrt(random_walk_sum[k]/L);
   		random_walk2[i][k] = random_walk[i][k]*random_walk[i][k];
   		
   		}
   	
      }
   
      
   double sum_prog = 0;
   double su2_prog = 0;
   double err_prog = 0;
   
   for (int k=0; k<100; k++) {
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog = 0;
   	su2_prog = 0;
   	err_prog = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog += random_walk[j][k];
   		su2_prog += random_walk2[j][k];
   	}
   	
   	sum_prog /= (i+1); // Cumulative average
    	su2_prog /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog = 0;
    		}
    	else {
    		err_prog = sqrt((su2_prog - sum_prog*sum_prog)/i); 
    		}
    	
    	output1 << k+1 << " " << L*i << " " << sum_prog << " " << err_prog << " " << endl;
    	
    	if (i==N-1) output2 << k+1 << " " << sum_prog << " " << err_prog << " " << endl;
	
	}
	
	}
	
	output1.close();
	output2.close();
	
	 // Esercizio 02.2.2: random walk in the continuum
	 
	ofstream output3;
   output3.open("random_walk_continuum_blocks.dat");
   
   ofstream output4;
   output4.open("random_walk_continuum.dat");

	double x1 = 0, x2 = 0; // Auxiliar variables to sample the unit sphere uniformly
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	fill(random_walk_sum.begin(), random_walk_sum.end(), 0);
   	
   	for (int j=0; j<L; j++) { // Iteration in each block
   		fill(position.begin(), position.end(), 0);
   		
   		for (int k=0; k<N_steps_rw; k++) {
   			
   			do {
   			x1 = rnd.Rannyu(-1,1);
   			x2 = rnd.Rannyu(-1,1);
   			} while (x1*x1+x2*x2>1);
   			 			
   			position[0] += a*2*x1*sqrt(1-x1*x1-x2*x2);
   			position[1] += a*2*x2*sqrt(1-x1*x1-x2*x2);
   			position[2] += a*(1-2*(x1*x1+x2*x2));
   			
   			random_walk_sum[k] += position[0]*position[0]+position[1]*position[1]+position[2]*position[2]; 
   		}
   	}
   	
   	for (int k=0; k<100; k++) {
   	
   		random_walk[i][k] = sqrt(random_walk_sum[k]/L);
   		random_walk2[i][k] = random_walk[i][k]*random_walk[i][k];
   		
   		}
   	
      }
      
   for (int k=0; k<100; k++) {
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
		sum_prog = 0;
   	su2_prog = 0;
   	err_prog = 0;
   	
   	for (int j=0; j<i+1;  j++) {
   		sum_prog += random_walk[j][k];
   		su2_prog += random_walk2[j][k];
   	}
   	
   	sum_prog /= (i+1); // Cumulative average
    	su2_prog /= (i+1); // Cumulative square average
    	
    	if (i == 0) {
    		err_prog = 0;
    		}
    	else {
    		err_prog = sqrt((su2_prog - sum_prog*sum_prog)/i); 
    		}
    	
    	output3 << k+1 << " " << L*i << " " << sum_prog << " " << err_prog << " " << endl;
    	
    	if (i==N-1) output4 << k+1 << " " << sum_prog << " " << err_prog << " " << endl;
	
	}
	
	}
	
	output3.close();
	output4.close();
   
   
	return 0;
}

