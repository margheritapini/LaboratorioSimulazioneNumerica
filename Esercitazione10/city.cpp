#include <iostream>
#include <math.h>
#include "city.h"

using namespace std;

City :: City() {
	_x.set_size(2);
}

int City :: getindex() {
	return _index;
}

void City :: setindex(int index) {
	_index = index;
}

double City :: getposition(int dim) {
	return _x(dim);
}

void City :: setposition(int dim, double position) {
	_x(dim) = position;
}

