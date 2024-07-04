#ifndef __funzioni_h__
#define __funzioni_h__

#include <cmath>
#include <vector>

using namespace std;

// Classe astratta FunzioneBase
class FunzioneBase {

    public:
        virtual double Eval(const vector<double>& x) const = 0;

        virtual ~FunzioneBase() {};
};

// Classe derivata concreta: Funzione d'onda atomo di idrogeno 100
class hydrogenatom_100: public FunzioneBase {

    public:

        // Costruttori
        hydrogenatom_100() {};
        // Distruttore
        virtual ~hydrogenatom_100() {};

        // Metodo per valutare la funzione in x
        virtual double Eval(const vector<double>& x) const {
            return 1.0 / M_PI * exp(-2*sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
        }

};

// Classe derivata concreta: Funzione d'onda atomo di idrogeno 210
class hydrogenatom_210: public FunzioneBase {

    public:

        // Costruttori
        hydrogenatom_210() {};
        // Distruttore
        virtual ~hydrogenatom_210() {};

        // Metodo per valutare la funzione in x
        virtual double Eval(const vector<double>& x) const {
            return 1.0 / 64.0 * 2.0 / M_PI * x[2] * x[2] * exp(- sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
        }

};

#endif

