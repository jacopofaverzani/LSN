#include <cmath>
#include <fstream>
#include <iostream>
#include "random.h"
using namespace std;

// Classe astratta per definire una funzione generica
class FunzioneBase_3D {
public:
    virtual double Val(double x[3]) const = 0; // Metodo virtuale puro per valutare la funzione in un punto
};

// Classe che rappresenta la densità di probabilità dello stato fondamentale
class stato_fondamentale: public FunzioneBase_3D {
public:
	stato_fondamentale() {;} 
	~stato_fondamentale() {;} 
	virtual double Val(double posizione[3]) const override {
        double r = sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
		return pow((1/(sqrt(M_PI)))*exp(-r), 2);
	};
};

// Classe che rappresenta la densità di probabilità dello stato eccitato 2p
class stato_eccitato: public FunzioneBase_3D {
public:
	stato_eccitato() {;} 
	~stato_eccitato() {;} 
	virtual double Val(double posizione[3]) const override {
        double r = sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
		return pow((1./8)*sqrt(2/M_PI)*exp(-r/2)*posizione[2], 2);
	};
};

// Algoritmo di Metropolis con probabilità di transizione uniforme
unsigned int Metropolis_unif(double posizione[3], double dim_step, Random& rnd, const FunzioneBase_3D *p, 
                            unsigned int passi_accettati);

// Algoritmo di Metropolis con probabilità di transizione normale multivariata
unsigned int Metropolis_multi(double posizione[3], double dim_step, Random& rnd, const FunzioneBase_3D *p, 
                            unsigned int passi_accettati);

// Funzione per calcolare l'incertezza statistica
double errore(double media, double media2, int n);
// Prende in ingresso media, media quadratica e numero di blocchi usati (n = N-1)