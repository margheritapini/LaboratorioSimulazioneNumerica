#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include "functions.h"
#include "Metropolis.h"

using namespace std;


int main (int argc, char *argv[]) {

    // Open output file to save the Hamiltonian values during simulated annealing
	ofstream output_H;
	output_H.open("H_SA.dat");

    // Instantiate Metropolis object and initialize it
	Metropolis Metr;
    Metr.initialize();
    
    // Define the initial and final values of beta (inverse temperature) and the increment for beta at each SA step
    double beta_in = 1; // T = 1 initial temperature
    double beta_fin = 500; // T = 1/beta_fin final temperature
    double beta = beta_in; // To instantly store the temperature 
    double beta_increase = 1;
    
    // Calculate the number of simulated annealing steps
    int SA_steps = static_cast<int>((beta_fin - beta_in)/beta_increase); // Simulated annealing steps
    
    int n_i; // Number of iterations for each SA step

    // PARAMETERS 
    double mu = 1; 
    double sigma = 1; 
  
    // Variables to store the mean and error of the Hamiltonian
    double mean_H; // To store the instant expectation value of the Hamiltonian on the current trial state
    double error_mean_H; // Statistical uncertainty associated to it
    // Define the number of blocks and throws for computing mean_H through Metropolis algorithm at each change of mu and sigma in SA iterations
    int n_blocks_metr_H = 100;
    int n_throws_metr_H = 100000;
    
    /*
    // Initial state: Compute the initial value of mean_H 
    Metr.compute_expect_H(mu, sigma, n_blocks_metr_H, n_throws_metr_H, mean_H, error_mean_H, false);
  
    // Loop over the simulated annealing steps
    for (int i = 0; i < SA_steps; i++) {
    	
        if (i != 0) beta += beta_increase; // Increase beta after the first iteration
    	
        // Calculate the number of iterations for this SA step
        n_i = 100 * pow(beta, 0.6);
            
        // Perform simulated annealing step
        Metr.SA(beta, mu, sigma, mean_H, error_mean_H, n_i, n_blocks_metr_H, n_throws_metr_H);
    	
        // Output the current state (current mean_H, error_mean_H, and mu and sigma which give those values) to file
        output_H << beta << " " << setprecision(18) << mean_H << " " << setprecision(18) << error_mean_H << " " << mu << " " << sigma << endl;
    }

    // Close the output file
	output_H.close();
	*/
	 // Optimized mu and sigma
	mu = 0.80;
	sigma = 0.62;
   Metr.compute_expect_H(mu, sigma, n_blocks_metr_H, n_throws_metr_H, mean_H, error_mean_H, true);
    
    return 0;
}
       
