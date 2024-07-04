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
    Metropolis() {;}; // Constructor
    ~Metropolis(){;}; // Destructor
    
    // Methods
    void initialize();
    vector<double> Unif(FunzioneBase* p, vector<double>& x_n, double delta);
    vector<double> Gaussians(FunzioneBase* p, vector<double>& x_n, double sigma);
    double choose_delta_Unif(FunzioneBase* p, vector<double>& x_n, double delta_guess, double acceptance);
    double choose_sigma_Gaussians(FunzioneBase* p, vector<double>& x_n, double sigma_guess, double acceptance);
    

private: 
    double alpha;
    double r;
    vector<double> x_new;
    Random T;
};

void Metropolis::initialize() {

	x_new.resize(3);

	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	int seed[4]; // Read the seed of the RNG
	ifstream Seed("seed.in");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	T.SetRandom(seed,p1,p2);
  
	Seed.close();
  }

vector<double> Metropolis::Unif(FunzioneBase* p, vector<double>& x_n, double delta) {
    x_new[0] = x_n[0] + T.Rannyu(-delta, delta);
    x_new[1] = x_n[1] + T.Rannyu(-delta, delta);
    x_new[2] = x_n[2] + T.Rannyu(-delta, delta);
    
    alpha = min(1.0, p->Eval(x_new) / p->Eval(x_n));
    
    r = T.Rannyu();
    
    if (alpha >= r) 
        return x_new;
    else 
        return x_n;
}

vector<double> Metropolis::Gaussians(FunzioneBase* p, vector<double>& x_n, double sigma) {
    x_new[0] = T.Gauss(x_n[0],sigma);
    x_new[1] = T.Gauss(x_n[1],sigma);
    x_new[2] = T.Gauss(x_n[2],sigma);
    
    alpha = min(1.0, p->Eval(x_new) / p->Eval(x_n));
    
    r = T.Rannyu();
    
    if (alpha >= r) 
        return x_new;
    else 
        return x_n;
}

double Metropolis::choose_delta_Unif(FunzioneBase* p, vector<double>& x_n, double delta_guess, double acceptance) {

    // Variables initialization
    int total_delta_trials_more = 0;
    int total_delta_trials_less = 0;
    int accepted = 0;
    int total_throws = 100;
    double delta;
    double resize = delta_guess / 50.0;
    
    if (delta_guess <= 0) {
    	cout << "The step for the uniform distribution must be larger than zero." << endl;
      return 0;
      }
    
    // Loop to find appropriate delta
    do {
        accepted = 0;

        // Adjust delta based on trial counts
        if ((total_delta_trials_more + total_delta_trials_less) % 2 == 0) {
            delta = delta_guess + resize * total_delta_trials_more;
            total_delta_trials_more++;
        } else {
            delta = delta_guess - resize * total_delta_trials_less;
            total_delta_trials_less++;
        }

        // Loop for throwing trials
        for (int i = 0; i < total_throws; i++) {
            // Generate new position
            x_new[0] = x_n[0] + T.Rannyu(-delta, delta);
            x_new[1] = x_n[1] + T.Rannyu(-delta, delta);
            x_new[2] = x_n[2] + T.Rannyu(-delta, delta);

            // Calculate acceptance probability
            alpha = min(1.0, p->Eval(x_new) / p->Eval(x_n));
            r = T.Rannyu();

            // Accept or reject move
            if (alpha >= r) {
                accepted++;
                x_n = x_new;
            }
        }

        // If the delta trials are exhausted, return 0
        if (total_delta_trials_less == 49) {
            cout << "Couldn't find a step for the uniform distribution with the requested acceptance. Try another initial guess." << endl;
            return 0;
        }

    // Loop until desired acceptance rate is achieved or delta trials are exhausted
    } while (double(accepted) / total_throws <= acceptance);
	
	cout << "accettanza: " << double(accepted) / total_throws << endl;
	
    return delta;
}

double Metropolis::choose_sigma_Gaussians(FunzioneBase* p, vector<double>& x_n, double sigma_guess, double acceptance) {

    // Variables initialization
    int total_sigma_trials_more = 0;
    int total_sigma_trials_less = 0;
    int accepted = 0;
    int total_throws = 100;
    double sigma;
    double resize = sigma_guess / 50.0;
    
    if (sigma_guess <= 0) {
    	cout << "The sigma of the Gaussian distribution must be larger than zero." << endl;
      return 0;
      }
    
    // Loop to find appropriate sigma
    do {
        accepted = 0;

        // Adjust sigma based on trial counts
        if ((total_sigma_trials_more + total_sigma_trials_less) % 2 == 0) {
            sigma = sigma_guess + resize * total_sigma_trials_more;
            total_sigma_trials_more++;
        } else {
            sigma = sigma_guess - resize * total_sigma_trials_less;
            total_sigma_trials_less++;
        }

        // Loop for throwing trials
        for (int i = 0; i < total_throws; i++) {
            // Generate new position
            x_new[0] = T.Gauss(x_n[0],sigma);
    			x_new[1] = T.Gauss(x_n[1],sigma);
    			x_new[2] = T.Gauss(x_n[2],sigma);

            // Calculate acceptance probability
            alpha = min(1.0, p->Eval(x_new) / p->Eval(x_n));
            r = T.Rannyu();

            // Accept or reject move
            if (alpha >= r) {
                accepted++;
                x_n = x_new;
            }
        }

        // If the sigma trials are exhausted, return 0
        if (total_sigma_trials_less == 49) {
            cout << "Couldn't find a sigma for the gaussian distribution with the requested acceptance. Try another initial guess." << endl;
            return 0;
        }

    // Loop until desired acceptance rate is achieved or sigma trials are exhausted
    } while (double(accepted) / total_throws <= acceptance);
		
		cout << "accettanza: " << double(accepted) / total_throws << endl;
    return sigma;
}


#endif


