#include <cmath>
#include <cstdlib>
#include <string>
#include "Genetic_Algorithm.h"

using namespace std;
using namespace arma;

void GA :: set_ncity(int ncity) {
	_ncity = ncity;
	}
	
int GA :: get_ncity() {
	return _ncity;
}
	
void GA :: set_npopulation(int npopulation) {
	_npopulation = npopulation;
	}
	
int GA :: get_npopulation() {
	return _npopulation;
}
	
void GA :: initialize(int rank) {

	int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("Primes");
  
  for (int i = 0; i <= rank; i++) {
        Primes >> p1 >> p2;
    }
    
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);
	
	ifstream input("input.dat"); // Start reading 
  string property;
  
  while ( !input.eof() ){
    input >> property;
    if( property == "N_CITIES" ){
      input >> _ncity;      
    } else if( property == "N_POPULATION" ){
      input >> _npopulation;
    } else if( property == "ENDINPUT" ){
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  input.close();
  
  _cities.set_size(_ncity);
  _population.set_size(_npopulation);
  
  for (int i = 0; i < _npopulation; i++) {
            _population(i) = individual(_ncity); // Initialize individual with the number of cities
        }
	}
	
void GA::initialize_positions() {
    ifstream input_cities("cap_prov_ita.dat");
    
    if (!input_cities.is_open()) {
        cerr << "Error opening file: cap_prov_ita.dat" << endl;
        return;
    }
    
    double x, y;
    
    input_cities >> fixed >> setprecision(10);
    
    for (int i = 0; i < _ncity; i++) {
        input_cities >> x >> y;
        _cities(i).setposition(0, x);
        _cities(i).setposition(1, y);
        _cities(i).setindex(i + 1);
    }

    input_cities.close();
}

void GA :: set_population() {

	for (int i=0; i<_ncity; i++) {
		_population(0).set_cityorder (i, i+1);
		}

	_population(0).measure_distance(_cities);
	
	for (int i=1; i<_npopulation; i++) {
		_population(i).permutation_cityorder(_population(i-1), int(_ncity/2), _rnd);
		_population(i).measure_distance(_cities);
	}
}

void GA :: print_population() const {

   /* for (int i = 0; i < _population.n_elem; ++i) {
        cout << "Individual " << i << ":" << endl;
        cout << "City order: " << _population(i)._chromosome.t(); // Transpose to print as row vector
        cout << "Distance: " << _population(i)._distance << endl;
    }*/
    
    cout << "Individual " << 0 << ": " << endl;
    cout << "City order: " << _population(0)._chromosome.t(); // Transpose to print as row vector
    cout << "Distance: " << _population(0)._distance << endl;
   
}

void GA :: save_path(int rank) {
	
	ofstream output;
	output.open("best_path" + to_string(rank) + ".dat");
	
	for (int i=0; i<_ncity; i++) {
		output << (_population(0).get_chromosome())(i) << endl;	
	}
	
	output << endl;
	output << _population(0)._distance;
	output.close();
}

void GA :: sort_population() {
	field <individual> sorted_population(_npopulation);
	vector <int> indices(_npopulation);
	
	for (int i=0; i<_npopulation; i++) {
		indices[i] = i;
	}
	
	sort(indices.begin(), indices.end(), [this](int i1, int i2) {
		return _population(i1)._distance < _population(i2)._distance; });
		
	for (int i=0; i<_npopulation; i++) {
		sorted_population(i) = _population(indices[i]);
		}
		
	_population = sorted_population;

}

void GA::initialize_output_file(int rank) {
    ofstream output("Best_distance" + to_string(rank) + ".dat");
    output.close();
}

void GA :: measure(int generation, int rank) {
	
	vec best_distances(static_cast<int>(_npopulation/2.));
	
	for(int i=0; i<static_cast<int>(_npopulation/2.); i++) {
		best_distances(i) = _population(i).get_distance();
		}
	
	ofstream output;
	output.open("Best_distance" + to_string(rank) + ".dat", ios::app);
	output << generation << setw(25) << _population(0)._distance << setw(18) << mean(best_distances) << endl;
	output.close();
	
}

int GA :: selection_operator() {
	double r = abs(_rnd.Rannyu());
	
	return int((_npopulation-1)/pow(0.5,8)*pow(r-0.5,8));
}

individual GA::_crossover(ivec genitore1, ivec genitore2, int quale_figlio) {
    
    ivec figlio(_ncity);
    
    ivec genitore1_copia = genitore1;
    ivec genitore2_copia = genitore2;
    
    int num;
    
    // int cut = static_cast<int>(2./5.* _ncity);
    int cut = static_cast<int>(abs(_rnd.Rannyu(1./5.* _ncity, 4./5.*_ncity)));
    
    if (quale_figlio == 0) {
        for (int i = 0; i < cut; i++) {
            num = genitore1(i);
            figlio(i) = num;
            genitore2_copia.replace(num, 0);
        }
        // Rimuovi gli zeri da genitore2_copia
        uvec non_zero_indices = find(genitore2_copia != 0);
        genitore2_copia = genitore2_copia(non_zero_indices);

        for (int i = cut; i < _ncity; i++) {
            figlio(i) = genitore2_copia(i - cut);
        }
    } else {
        for (int i = 0; i < cut; i++) {
            num = genitore2(i);
            figlio(i) = num;
            genitore1_copia.replace(num, 0);
        }
        // Rimuovi gli zeri da genitore1_copia
        uvec non_zero_indices = find(genitore1_copia != 0);
        genitore1_copia = genitore1_copia(non_zero_indices);
        for (int i = cut; i < _ncity; i++) {
            figlio(i) = genitore1_copia(i - cut);
        }
    }
    
    // Creazione dell'individuo e calcolo della distanza misurata
    individual i_figlio(_ncity, figlio);
    
    return i_figlio;
}

individual GA :: _mutation(individual figlio) {

	double m1 = _rnd.Rannyu();
	double m2 = _rnd.Rannyu();
	double m3 = _rnd.Rannyu();
	double m4 = _rnd.Rannyu();
	double m5 = _rnd.Rannyu();	
	
	// PERMUTATION OF TWO RANDOM CITIES
	if (m1 < 0.021) figlio.permutation_cityorder(figlio, 1, _rnd);
	
	
	// PERMUTATION OF M CONTIGUOUS CITIES WITH OTHER M CONTIGUOUS CITIES
	if (m2 < 0.21) {
		
		int block_size = static_cast<int>(_rnd.Rannyu(1, _ncity/2.-1));
		int start1_index = static_cast<int>(_rnd.Rannyu(1, _ncity - 2*block_size));
		int start2_index = static_cast<int>(_rnd.Rannyu(start1_index + block_size, _ncity - block_size));

		int temp;
		
		for (int i = 0; i < block_size; i++) {
		
			temp = (figlio.get_chromosome())(start1_index+i);
			figlio.set_gene(start1_index+i, (figlio.get_chromosome())(start2_index+i));
			figlio.set_gene(start2_index+i, temp);
		}
		} 
	
	// INVERSION OF A RANDOMIC SUB-BLOCKS OF CITIES
	
	if (m3 < 0.191) {
    int block_size = static_cast<int>(_rnd.Rannyu(1, _ncity-1));
    int start_index = static_cast<int>(_rnd.Rannyu(1, _ncity - block_size));
    int end_index = start_index + block_size - 1;

    for (int i = 0; i < block_size / 2; i++) {
        int temp = (figlio.get_chromosome())(start_index + i);
        figlio.set_gene(start_index + i, (figlio.get_chromosome())(end_index - i));
        figlio.set_gene(end_index - i, temp);
    }
}
	
	// N (RANDOMIC) POSITIONS - SHIFT FOR ALL CITIES 

	if (m4 < 0.161) {
		ivec copy(_ncity-1);
		ivec shifted_copy(_ncity-1);
		
		for (int i=0; i<_ncity-1; i++) {
			copy(i) = (figlio.get_chromosome())(i+1);}
			
		shifted_copy = shift(copy, _rnd.Rannyu(1, _ncity-1));
		
		for (int i=0; i<_ncity-1; i++) {
			figlio.set_gene(i+1, shifted_copy(i));
			}
	}
	
	// N (RANDOMIC) POSITIONS - SHIFT FOR A SUB-GROUP OF CITIES 
	if (m5 < 0.11) {
		
		int block_size = static_cast<int>(_rnd.Rannyu(1, _ncity/2.-1));
		int start1_index = static_cast<int>(_rnd.Rannyu(1, _ncity - 2*block_size));
		int start2_index = static_cast<int>(start1_index + block_size);

		int temp;
		
		for (int i = 0; i < block_size; i++) {
		
			temp = (figlio.get_chromosome())(start1_index+i);
			figlio.set_gene(start1_index+i, (figlio.get_chromosome())(start2_index+i));
			figlio.set_gene(start2_index+i, temp);
		}
		} 
	
	return figlio;
}

void GA :: new_generation() {

	field <individual> new_population = _population;
	
	double p_c = 0.85; // 1
	
	for (int i=0; i<_npopulation; i+=2) {
		int j1 = selection_operator();
		int j2 = selection_operator();
		
		individual genitore1 = _population(j1);
		individual genitore2 = _population(j2);
		
		if (_rnd.Rannyu() < p_c) new_population(i) = _crossover(genitore1.get_chromosome(), genitore2.get_chromosome(), 0);
		if (_rnd.Rannyu() < p_c) new_population(i+1) = _crossover(genitore1.get_chromosome(), genitore2.get_chromosome(), 1);
		
		new_population(i) = _mutation(new_population(i));
		new_population(i+1) = _mutation(new_population(i+1));
		
		new_population(i).measure_distance(_cities);
		new_population(i+1).measure_distance(_cities);
		
		}
		
	_population = new_population;
	
	}
	
	

vector <int> GA :: choose_couples (int size) {
	vector <int> couples(size);
	iota(couples.begin(), couples.end(), 0);
	
	int j1, j2;
	int aux;
		
	for (int i=0; i<int(size/2); i++) {
		j1 = abs(_rnd.Rannyu(0, size));
		j2 = abs(_rnd.Rannyu(0, size));
			
			
		aux =couples[j1];
		couples[j1] =couples[j2];
		couples[j2] = aux;
		}
		
	for (int i=0; i<size; i++) {
		if (couples[i] == i) {
			if (i != (size-1)) {
				aux = couples[i];
				couples[i] = couples[i+1];
				couples[i+1] = aux;
				}
			else {
				aux = couples[i];
				couples[i] = couples[i-1];
				couples[i-1] = aux;
			}
			}
		}
		
	return couples;
}

int* GA :: get_best_chromosome() {
    
    int* best_chromosome = new int[_ncity];

    for (int i = 0; i < _ncity; i++) {
        best_chromosome[i] = (_population(0).get_chromosome())(i);
    }

    return best_chromosome;
}

void GA :: set_best_chromosome(int* best_chromosome) {

	for (int i=0; i<_ncity; i++) {
		_population(0).set_gene(i, best_chromosome[i]);
	}
	
	_population(0).measure_distance(_cities);

}

