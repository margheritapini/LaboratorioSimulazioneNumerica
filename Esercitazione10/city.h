#ifndef __City__
#define __City__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class City {

private:
  int _index;           // Index of the city
  vec _x;               // Position of the city vector

public: // Function declarations
  City();
  int  getindex();                        // Get the spin of the particle
  void setindex(int index);                // Set the spin of the particle
  double getposition(int dim);// Get the position of the particle along a specific dimension
  void   setposition(int dim, double position); // Set the position of the particle along a specific dimension
};

#endif // __City__
