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
#include "functions.h"
#include "Metropolis.h"

using namespace std;

double error(double acc, double acc2, int blk) {
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int main (int argc, char *argv[]) {

	 Metropolis Metr;
    Metr.initialize();
    
    // Vectors for hydrogen_100 and hydrogen_210
    vector<double> x_n_100(3, 0);
    vector<double> x_n_210(3, 0);

    FunzioneBase *hydrogen_100 = new hydrogenatom_100();
    FunzioneBase *hydrogen_210 = new hydrogenatom_210();
    
    // Fixing initial position for hydrogen_100
    x_n_100[0] = 0;
    x_n_100[1] = 1;
    x_n_100[2] = 1;

    // Fixing initial position for hydrogen_210
    x_n_210[0] = 1;
    x_n_210[1] = 1;
    x_n_210[2] = 5;
    
    int M = 2000000; // Total number of throws
    int N = 200;      // Number of blocks
    int L = int(M/N); // Number of throws in each block
    
    int M_eq = 10000;
    
    double block_av_100 = 0;
    double block_av_210 = 0;
    double average_100 = 0;
    double average_210 = 0;
    double global_av_100 = 0;
    double global_av_210 = 0;
    double global_av2_100 = 0;
    double global_av2_210 = 0;


	 // UNIFORM TRANSITION PROBABILITY
    ofstream cout_unif_r100;
    cout_unif_r100.open("OUTPUT/unif_r100.dat");

    ofstream cout_unif_coordinates_100;
    cout_unif_coordinates_100.open("OUTPUT/unif_coordinates_100.dat");

    ofstream cout_unif_r210;
    cout_unif_r210.open("OUTPUT/unif_r210.dat");

    ofstream cout_unif_coordinates_210;
    cout_unif_coordinates_210.open("OUTPUT/unif_coordinates_210.dat");

    double delta_guess_100 = 1.5;
    double delta_guess_210 = 3.2;
    double delta_100 = 0;
    double delta_210 = 0;

    // Finding delta for hydrogen_100
    delta_100 = Metr.choose_delta_Unif(hydrogen_100, x_n_100, delta_guess_100, 0.5);
    cout << "delta_100: " << delta_100 << endl;

    // Finding delta for hydrogen_210
    delta_210 = Metr.choose_delta_Unif(hydrogen_210, x_n_210, delta_guess_210, 0.5);
    cout << "delta_210: " << delta_210 << endl;
    
    for (int i=0; i<M_eq; i++) {
    	x_n_100 = Metr.Unif(hydrogen_100, x_n_100, delta_100);
    	x_n_210 = Metr.Unif(hydrogen_210, x_n_210, delta_210);
    	}
    
    // Iteration over the blocks for hydrogen_100 
    for (int i = 0; i < N; i++) {
        block_av_100 = 0;

        for (int j = 0; j < L; j++) {
            x_n_100 = Metr.Unif(hydrogen_100, x_n_100, delta_100);
            block_av_100 += sqrt(x_n_100[0]*x_n_100[0] + x_n_100[1]*x_n_100[1] + x_n_100[2]*x_n_100[2]);

            if (j % 100 == 0) {
                cout_unif_coordinates_100 << x_n_100[0] << " " << x_n_100[1] << " " << x_n_100[2] << endl;
            }
        }

        average_100 = block_av_100 / L;
        global_av_100 += average_100;
        global_av2_100 += average_100 * average_100;

        cout_unif_r100 << i + 1 << " " << global_av_100 / (i + 1) << " " << error(global_av_100, global_av2_100, i + 1) << endl;
    }

    // Iteration over the blocks for hydrogen_210
    for (int i = 0; i < N; i++) {
        block_av_210 = 0;

        for (int j = 0; j < L; j++) {
            x_n_210 = Metr.Unif(hydrogen_210, x_n_210, delta_210);
            block_av_210 += sqrt(x_n_210[0]*x_n_210[0] + x_n_210[1]*x_n_210[1] + x_n_210[2]*x_n_210[2]);

            if (j % 100 == 0) {
                cout_unif_coordinates_210 << x_n_210[0] << " " << x_n_210[1] << " " << x_n_210[2] << endl;
            }
        }

        average_210 = block_av_210 / L;
        global_av_210 += average_210;
        global_av2_210 += average_210 * average_210;

        cout_unif_r210 << i + 1 << " " << global_av_210 / (i + 1) << " " << error(global_av_210, global_av2_210, i + 1) << endl;
    }

    cout_unif_r100.close();
    cout_unif_coordinates_100.close();
    cout_unif_r210.close();
    cout_unif_coordinates_210.close();
    
    // GAUSSIAN TRANSITION PROBABILITY
    
    block_av_100 = 0;
    block_av_210 = 0;
    average_100 = 0;
    average_210 = 0;
    global_av_100 = 0;
    global_av_210 = 0;
    global_av2_100 = 0;
    global_av2_210 = 0;
    
    ofstream cout_gauss_r100;
    cout_gauss_r100.open("OUTPUT/gauss_r100.dat");

    ofstream cout_gauss_coordinates_100;
    cout_gauss_coordinates_100.open("OUTPUT/gauss_coordinates_100.dat");

    ofstream cout_gauss_r210;
    cout_gauss_r210.open("OUTPUT/gauss_r210.dat");

    ofstream cout_gauss_coordinates_210;
    cout_gauss_coordinates_210.open("OUTPUT/gauss_coordinates_210.dat");

    double sigma_guess_100 = 1;
    double sigma_guess_210 = 2;
    double sigma_100 = 0;
    double sigma_210 = 0;
    
    // Fixing initial position for hydrogen_100
    x_n_100[0] = 0;
    x_n_100[1] = 1;
    x_n_100[2] = 1;

    // Fixing initial position for hydrogen_210
    x_n_210[0] = 1;
    x_n_210[1] = 1;
    x_n_210[2] = 5;
    
    // Finding sigma for hydrogen_100
    sigma_100 = Metr.choose_sigma_Gaussians(hydrogen_100, x_n_100, sigma_guess_100, 0.5);
    
    cout << "sigma_100: " << sigma_100 << endl;

    // Finding sigma for hydrogen_210
    sigma_210 = Metr.choose_sigma_Gaussians(hydrogen_210, x_n_210, sigma_guess_210, 0.5);
    cout << "sigma_210: " << sigma_210 << endl;
    
    for (int i=0; i<M_eq; i++) {
    	x_n_100 = Metr.Gaussians(hydrogen_100, x_n_100, sigma_100);
    	x_n_210 = Metr.Gaussians(hydrogen_210, x_n_210, sigma_210);
    	}

    // Iteration over the blocks for hydrogen_100
    for (int i = 0; i < N; i++) {
        block_av_100 = 0;

        for (int j = 0; j < L; j++) {
            x_n_100 = Metr.Gaussians(hydrogen_100, x_n_100, sigma_100);
            block_av_100 += sqrt(x_n_100[0]*x_n_100[0] + x_n_100[1]*x_n_100[1] + x_n_100[2]*x_n_100[2]);

            if (j % 100 == 0) {
                cout_gauss_coordinates_100 << x_n_100[0] << " " << x_n_100[1] << " " << x_n_100[2] << endl;
            }
        }

        average_100 = block_av_100 / L;
        global_av_100 += average_100;
        global_av2_100 += average_100 * average_100;

        cout_gauss_r100 << i + 1 << " " << global_av_100 / (i + 1) << " " << error(global_av_100, global_av2_100, i + 1) << endl;
    }

    // Iteration over the blocks for hydrogen_210
    for (int i = 0; i < N; i++) {
        block_av_210 = 0;

        for (int j = 0; j < L; j++) {
            x_n_210 = Metr.Gaussians(hydrogen_210, x_n_210, sigma_210);
            block_av_210 += sqrt(x_n_210[0]*x_n_210[0] + x_n_210[1]*x_n_210[1] + x_n_210[2]*x_n_210[2]);

            if (j % 100 == 0) {
                cout_gauss_coordinates_210 << x_n_210[0] << " " << x_n_210[1] << " " << x_n_210[2] << endl;
            }
        }

        average_210 = block_av_210 / L;
        global_av_210 += average_210;
        global_av2_210 += average_210 * average_210;

        cout_gauss_r210 << i + 1 << " " << global_av_210 / (i + 1) << " " << error(global_av_210, global_av2_210, i + 1) << endl;
    }

    cout_gauss_r100.close();
    cout_gauss_coordinates_100.close();
    cout_gauss_r210.close();
    cout_gauss_coordinates_210.close();
    
    return 0;
}
