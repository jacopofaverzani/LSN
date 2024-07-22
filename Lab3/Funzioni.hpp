#ifndef __Funzioni_hpp__
#define __Funzioni_hpp__

#include <cmath>
#include <vector>
#include "random.h" // Includi la libreria per i numeri casuali

using namespace std;

// Funzione per calcolare l'incertezza statistica
double err(double mean, double mean2, unsigned int n);
// Prende in ingresso media, media quadratica e numero di blocchi usati (n = N-1)

// Funzione per il campionamento diretto del prezzo finale dell'asset
double European_dir(double risk, double vol, double S0, double T, Random& rnd);
// Prende in ingresso fattore di rischio, volatilità, prezzo iniziale, tempo di scadenza, prezzo di esercizio e 
// il generatore di numeri casuali

// Funzione per il campionamento tramite il percorso discreto del prezzo finale dell'asset
double European_step(double risk, double vol, double S0, double T, Random& rnd);
// Prende in ingresso fattore di rischio, volatilità, prezzo iniziale, tempo di scadenza, prezzo di esercizio e 
// il generatore di numeri casuali

#endif // Fine del file delle intestazioni