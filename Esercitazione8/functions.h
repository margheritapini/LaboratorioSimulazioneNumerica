#ifndef __funzioni_h__
#define __funzioni_h__

#include <cmath>
#include <vector>

using namespace std;

// Abstract base class for functions
class FunzioneBase {
    public:
        // Pure virtual function to evaluate the function at x
        virtual double Eval(const double x, double mu, double sigma) const = 0;
        // Virtual destructor for proper cleanup of derived classes
        virtual ~FunzioneBase() {};
};

// Derived class for a specific trial wavefunction squared
class psi_trial2: public FunzioneBase {
    public:
        // Constructor
        psi_trial2() {};
        // Destructor
        virtual ~psi_trial2() {};
        // Method to evaluate the function at x
        virtual double Eval(const double x, double mu, double sigma) const {
            double sigma2 = sigma * sigma;
            double minus2 = (x - mu) * (x - mu);
            double plus2 = (x + mu) * (x + mu);
            double psi = exp(-minus2 / (2 * sigma2)) + exp(-plus2 / (2 * sigma2));
            return psi * psi; // Return the square of the trial wavefunction
        }
};


// Derived class for the Hamiltonian function
class H_x: public FunzioneBase {
    public:
        // Constructor
        H_x() {};
        // Destructor
        virtual ~H_x() {};
        // Method to evaluate the Hamiltonian at x
        virtual double Eval(const double x, double mu, double sigma) const {
            // Calculate the potential V(x)
            double V = (x * x - 2.5) * x * x;

            // Calculate the wavefunction terms psiplus and psiminus
            double psiplus = exp(-0.5 * pow((x + mu) / sigma, 2.0));
            double psiminus = exp(-0.5 * pow((x - mu) / sigma, 2.0));

            // Calculate the kinetic energy term T
            double term1 = pow((x - mu) / sigma, 2.0);
            double term2 = pow((x + mu) / sigma, 2.0);
            double T = -0.5 * (psiminus * (term1 - 1) / (sigma * sigma) + psiplus * (term2 - 1) / (sigma * sigma));

            // Return the Hamiltonian value
            return V + T / (psiminus + psiplus);
        }
};

#endif

