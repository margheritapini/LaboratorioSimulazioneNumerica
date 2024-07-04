#ifndef __GeneticAlgorithm__
#define __GeneticAlgorithm__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h>
#include "city.h"
#include "random.h"

using namespace std;
using namespace arma;

struct individual {
	
	int _ncity;
	double _distance;
	ivec _chromosome;	// Vettore di interi che tiene conto dell'ordine delle città
	
	individual() {}
	
	individual(int ncity) {
		_ncity = ncity;
		_chromosome.set_size(_ncity);
		_chromosome(0) = 1;
	}
	
	individual(int ncity, ivec chromosome) {
		_ncity = ncity;
		_chromosome.set_size(_ncity);
		_chromosome = chromosome;
	}
	
	void set_cityorder (int index, int index_city) {
		_chromosome(index) = index_city;
	}
	
	ivec get_chromosome () {
		return _chromosome;
	}
	
	void set_chromosome (ivec chromosome) {
		_chromosome = chromosome;
	}
	
	void set_gene (int index, int new_value) {
		_chromosome(index) = new_value;
	}
	
	void permutation_cityorder (individual other, int n_permutations, Random* _rnd) {
		int j1, j2;
		double aux;
		
		_chromosome = other.get_chromosome();
		
		for (int i=0; i<n_permutations; i++) {
			j1 = _rnd->Rannyu(1, _ncity-0.0001);
			j2 = _rnd->Rannyu(1, _ncity-0.0001);
			
			aux = _chromosome(j1);
			_chromosome(j1) = _chromosome(j2);
			_chromosome(j2) = aux;
		}
	}
	
	void measure_distance(field<City> cities_pos) {
    _distance = 0;

    for (int i = 0; i < _ncity - 1; i++) {
        
        double city1_x = cities_pos(_chromosome(i)-1).getposition(0);
        double city1_y = cities_pos(_chromosome(i)-1).getposition(1);
        double city2_x = cities_pos(_chromosome(i + 1)-1).getposition(0);
        double city2_y = cities_pos(_chromosome(i + 1)-1).getposition(1);

        _distance += sqrt(pow(city1_x - city2_x, 2) + pow(city1_y - city2_y, 2));
    }

    double first_city_x = cities_pos(_chromosome(0)-1).getposition(0);
    double first_city_y = cities_pos(_chromosome(0)-1).getposition(1);
    double last_city_x = cities_pos(_chromosome(_ncity - 1)-1).getposition(0);
    double last_city_y = cities_pos(_chromosome(_ncity - 1)-1).getposition(1);

    _distance += sqrt(pow(first_city_x - last_city_x, 2) + pow(first_city_y - last_city_y, 2));
}

	double get_distance() {
	return _distance; }
};

class GA {

private:

  int _ncity;           // Number of cities
  int _npopulation;		  // Total number of chromosomes (rank of population)
  int _city_pos;				// How to generate the positions of the cities
  
  Random _rnd;          // Random number generator
  field <City> _cities; // Field of city objects representing the geography
  field <individual> _population; // Field of individuals: the different paths are registered and the distance
  
  vec _population_distances;
 

public: // Function declarations

	void set_ncity(int ncity);
	void set_npopulation(int npopulation);
	void initialize();
	void initialize_positions();
	void set_population();
	void print_population() const;
	void save_path();
	void sort_population();
	void measure(int generation);
	
	int selection_operator();
	individual _crossover(ivec, ivec, int);
	individual _mutation(individual);
	void new_generation();

};

#endif // __GeneticAlgorithm__

