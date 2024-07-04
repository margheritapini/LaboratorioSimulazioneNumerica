#ifndef __Metropolis_h__
#define __Metropolis_h__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "functions.h"

using namespace std;

class Metropolis {
public:
    // Constructor
    Metropolis() {;}
    // Destructor
    ~Metropolis(){;}

    // Methods
    void initialize(); // Initialize the random number generator
    double sample_psi2(FunzioneBase* p, double mu, double sigma, double x_n, double delta); // Perform a single Metropolis step with sample_psi2orm proposal distribution to sample the probability distribution |psi_T|^2
    double choose_delta_sample_psi2(FunzioneBase* p, double mu, double sigma, double x, double delta_guess, double acceptance); // Choose delta to achieve desired acceptance rate for sample_psi2 method
    bool sample_SA(double beta, double H_n, double H_new); // Decide whether to accept a new trial state based on Metropolis criterion
    void SA(double beta, double& mu, double& sigma, double& H, double& error_H, int n_i, int n_blocks, int n_throws); // Simulated annealing process
    void compute_expect_H(double mu, double sigma, int n_blocks, int n_throws, double& H, double& error_mean_H, bool final); // Compute expected value of Hamiltonian on a given psi_T state
    double error(double acc, double acc2, int blk); // Compute statistical error

private: 
    double alpha;
    Random T;
};

// Initialize the random number generator
void Metropolis::initialize() {
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    int seed[4]; // Read the seed of the RNG
    ifstream Seed("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    T.SetRandom(seed, p1, p2);
    Seed.close();
}

// Perform a single Metropolis step with sample_psi2orm proposal distribution to sample the probability distribution |psi_T|^2
double Metropolis::sample_psi2(FunzioneBase* p, double mu, double sigma, double x_n, double delta) {
    double x_new = x_n + T.Rannyu(-delta, delta);
    alpha = min(1.0, p->Eval(x_new, mu, sigma) / p->Eval(x_n, mu, sigma));

    if (alpha >= T.Rannyu()) 
        return x_new;
    else 
        return x_n;
}

// Decide whether to accept a new trial state based on Metropolis criterion
bool Metropolis::sample_SA(double beta, double H_n, double H_new) {
    double delta_E = H_new - H_n;
    alpha = min(1.0, exp(-beta * delta_E));

    if (alpha >= T.Rannyu()) 
        return true;
    else 
        return false;
}

// Choose delta to achieve desired acceptance rate for sample_psi2 method
double Metropolis::choose_delta_sample_psi2(FunzioneBase* p, double mu, double sigma, double x, double delta_guess, double acceptance) {
    int total_delta_trials_more = 0;
    int total_delta_trials_less = 0;
    int accepted = 0;
    int total_throws = 100;
    double delta;
    double resize = delta_guess / 50.0;
    
    double x_new;

    if (delta_guess <= 0) {
        cout << "The step for the sample_psi2 distribution must be larger than zero." << endl;
        return 0;
    }

    do {
        accepted = 0;
        if ((total_delta_trials_more + total_delta_trials_less) % 2 == 0) {
            delta = delta_guess + resize * total_delta_trials_more;
            total_delta_trials_more++;
        } else {
            delta = delta_guess - resize * total_delta_trials_less;
            total_delta_trials_less++;
        }

        for (int i = 0; i < total_throws; i++) {
            x_new = x + T.Rannyu(-delta, delta);
            alpha = min(1.0, p->Eval(x_new, mu, sigma) / p->Eval(x, mu, sigma));
            
            if (alpha >= T.Rannyu()) {
                accepted++;
                x = x_new;
            }
        }

        if (total_delta_trials_less == 10000) {
            cout << "Couldn't find a step for the sample_psi2 distribution with the requested acceptance. Delta_guess will be used." << endl;
            return delta_guess;
        }

    } while (!((acceptance - 0.05) <= (static_cast<double>(accepted) / static_cast<double>(total_throws)) && (static_cast<double>(accepted) / static_cast<double>(total_throws)) <= (acceptance + 0.05)));
	
    return delta;
}

// Simulated annealing process
void Metropolis::SA(double beta, double& mu, double& sigma, double& H, double& error_H, int n_i, int n_blocks, int n_throws) {
    double mu_prop, sigma_prop, H_prop, error_H_prop;

    for (int i = 0; i < n_i; i++) { // At each beta, iterate n_i times the mu and sigma research
    	  // Proposed mu and sigma
        mu_prop = abs(mu + pow(beta, -0.5) * T.Rannyu(-1, 1));
        sigma_prop = abs(sigma + pow(beta, -0.5) * T.Rannyu(-1, 1));
        // Compute <H> on the new trial state 
        compute_expect_H(mu_prop, sigma_prop, n_blocks, n_throws, H_prop, error_H_prop, false);
        // Check if the new state is accepted adn, if yes, change mu and sigma with the new ones
        if (sample_SA(beta, H, H_prop)) {
            mu = mu_prop;
            sigma = sigma_prop;
            H = H_prop;
            error_H = error_H_prop;
        }
    }
}

// Compute expected value of Hamiltonian
void Metropolis::compute_expect_H(double mu, double sigma, int n_blocks, int n_throws, double& H, double& error_mean_H, bool final) {
    ofstream output_final_H;
    ofstream output_final_coordinates;
    if (final) {
        output_final_H.open("H_final.dat");
        output_final_coordinates.open("coordinates.dat");
    }

    int M = n_throws; // Total number of throws
    int N = n_blocks;  // Number of blocks
    int L = int(M / N); // Number of throws in each block
    
    int M_eq = 1000; // Equilibration throws

    double block_av = 0;
    double average = 0;
    double global_av = 0;
    double global_av2 = 0;

    double x = 0; // Initial position
    double delta = 2;

    FunzioneBase *mean_H = new H_x();
    FunzioneBase *Tpsi2 = new psi_trial2();

    delta = choose_delta_sample_psi2(Tpsi2, mu, sigma, x, delta, 0.5);

    for (int i = 0; i < M_eq; i++) {
        x = sample_psi2(Tpsi2, mu, sigma, x, delta);
    }

    for (int i = 0; i < N; i++) {
        block_av = 0;

        for (int j = 0; j < L; j++) {
            x = sample_psi2(Tpsi2, mu, sigma, x, delta);
            block_av += mean_H->Eval(x, mu, sigma);
            
            if (final==true and j % 10 == 0) {
                output_final_coordinates << x << endl;
        }
        }

        average = block_av / L;
        global_av += average;
        global_av2 += average * average;

        if (final) {
            output_final_H << i + 1 << " " << global_av / (i + 1) << " " << error(global_av, global_av2, i + 1) << endl;
        }
    }

    if (final) {
        output_final_H.close();
        output_final_coordinates.close();
    }

    delete mean_H;
    delete Tpsi2;

    error_mean_H = error(global_av, global_av2, N);
    H = global_av / N;
}

// Compute statistical error
double Metropolis::error(double acc, double acc2, int blk) {
    if (blk <= 1) return 0.0;
    else return sqrt(fabs(acc2 / double(blk) - pow(acc / double(blk), 2)) / double(blk));
}

#endif

