#ifndef __funzioni__
#define __funzioni__

#include "random.h"

// Funzione per calcolare l'incertezza statistica
double err(double media, double media2, unsigned int n);

// Funzione per calcolare il chi quadrato
double chi2(unsigned int oss, int teo);

// Funzione per campionare un angolo tra 0 e pi
double angle(Random& rnd);

#endif