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
#include "city.h"
#include "Genetic_Algorithm.h"

using namespace std;

int main (int argc, char *argv[]) {

	GA GenAlg;
	GenAlg.initialize();
	GenAlg.initialize_positions();
	GenAlg.set_population();
	GenAlg.sort_population();
	
	for (int i=0; i<1000; i++) {
		GenAlg.new_generation();
		GenAlg.sort_population();
		GenAlg.measure(i);
		}

	GenAlg.print_population();
	GenAlg.save_path();
	
	return 0;
}
