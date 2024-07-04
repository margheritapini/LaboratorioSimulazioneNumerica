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
   
   double L = 9; // Lenght of the needle
   double d = 10; // Distance between lines
   
   double y = 0; // y position of one of the two endpoint of the needle
   double phi = 0; // Slope of the needle
   
   int n = 0; // int variable to express the position of y in terms of d: y = n*d + const (const<d), to identify the vertical d-lenght "bands" in which the plane is divided
   
   int counts = 0; // Number of times the needle hits the lines
   
   ofstream output;
   output.open("pi_values.dat");
   
   int M = 1000000; // Total number of throws
   int N = 100; // Number of blocks
   int block_throws = int(M/N); // Number of throws in each block
   
   double sum_prog = 0;
   double su2_prog = 0;
   double err_prog = 0;
   
   double pi;
   
   for (int i=0; i<N; i++) { // Iteration over the blocks
   	counts = 0;
   	
   	for (int j=0; j<block_throws; j++) { // Iteration in each block
   		
   		// Position of the needle: expressed in terms of one endpoint and the angle
   		y = rnd.Rannyu(0, RAND_MAX);
   		phi = rnd.Rannyu(0, RAND_MAX);
   	   	
   		if (sin(phi)<0) y += L*sin(phi); // If (sin(phi)<0), y corresponds to the position of the top endpoint. I translate y in order to take it as the position of the bottom endpoin of the needle
   	
   		n = static_cast<int>(y/d); // I evaluate in which "band"
   	
   		if (y+L*abs(sin(phi))>d*(n+1)) counts++; // Condition of intersection
   	}
   	
   	pi = 2*L*block_throws/(counts*d);
   	
   	sum_prog += pi;
   	su2_prog += pi*pi;
   	
   	if (i == 0) {
    		err_prog = 0;
    		}
    	else {
    		err_prog = sqrt((su2_prog/(i+1) - pow(sum_prog/(i+1),2)) /i); 
    		}
		
		output << block_throws*i << " " << sum_prog/(i+1) << " " << err_prog << " " << endl;
   	
      }
	
	output.close();
	rnd.SaveSeed();
   
	return 0;
}
