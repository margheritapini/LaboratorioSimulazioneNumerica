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

#include "mpi.h"

using namespace std;


int main(int argc, char* argv[]) {

    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (size % 2 != 0) {
        cout << "Size must be even!" << endl;
        MPI_Finalize();
        return 0;
    }

    GA GenAlg;
    GenAlg.initialize(rank);
    GenAlg.initialize_positions();
    
    GenAlg.initialize_output_file(rank);

    GenAlg.set_population();
    GenAlg.sort_population();
    
    vector<int> couples(size);

    for (int i = 0; i < 3000; i++) {
        
        GenAlg.new_generation();
        GenAlg.sort_population();
			
        if (i%20 == 0) {
        		int * best_chromosome_send = GenAlg.get_best_chromosome();
        		int * best_chromosome_recv = new int[GenAlg.get_ncity()];
            
            // Genero un vettore con le coppie ordinate e le comunico
            if (rank == 0) couples = GenAlg.choose_couples(size);
            MPI_Bcast(couples.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
            
            // Scambia best_individual_send tra i processi in couples
   	 		for (int j = 0; j < size; j += 2) {
        			int partner1 = couples[j];
        			int partner2 = couples[j + 1];

        			if (rank == partner1) {
            MPI_Send(best_chromosome_send, GenAlg.get_ncity(), MPI_INT, partner2, 0, MPI_COMM_WORLD);
            MPI_Recv(best_chromosome_recv, GenAlg.get_ncity(), MPI_INT, partner2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else if (rank == partner2) {
            MPI_Recv(best_chromosome_recv, GenAlg.get_ncity(), MPI_INT, partner1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(best_chromosome_send, GenAlg.get_ncity(), MPI_INT, partner1, 0, MPI_COMM_WORLD);
        }
    	}
            
            // Aggiorna la popolazione con il miglior individuo ricevuto
            GenAlg.set_best_chromosome(best_chromosome_recv);
            GenAlg.sort_population();
        }
 
        GenAlg.measure(i, rank); 
    }

    GenAlg.print_population();
    GenAlg.save_path(rank);

    MPI_Finalize();

    return 0;
}
