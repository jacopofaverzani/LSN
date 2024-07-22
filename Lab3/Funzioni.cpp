#include "Funzioni.hpp"
	
// Funzione per calcolare l'incertezza statistica
double err(double mean, double mean2, unsigned int n){ 
    // Al primo blocco pone l'incertezza uguale a zero
    if(n == 0){
        return 0;
    } else { // Altrimenti, calcola e restituisce l'incertezza statistica
        return sqrt((mean2 - mean*mean) / n);
    }
}

// Funzione per il campionamento diretto del prezzo finale dell'asset
double European_dir(double risk, double vol, double S0, double T, Random& rnd){
    double z = rnd.Gauss(0,1);
    return S0*exp((risk - 0.5*vol*vol)*T + vol*z*sqrt(T));
}

// Funzione per il campionamento tramite il percorso discreto del prezzo finale dell'asset
double European_step(double risk, double vol, double S0, double T, Random& rnd){
    double S_t = S0;
    double delta_t = T/100;
    for(unsigned int i = 1; i < 100; i++){
        double z = rnd.Gauss(0,1);
        S_t *= exp((risk - 0.5*vol*vol)*delta_t + vol*z*sqrt(delta_t));
    }
    return S_t;
}

