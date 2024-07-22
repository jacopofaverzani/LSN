#include "funzioni.hpp" 
#include <cmath>

// Funzione per calcolare l'incertezza statistica
double err(double media, double media2, unsigned int n){ 
    if(n == 0) return 0;                               // Al primo blocco pone l'incertezza uguale a zero
    else return sqrt((media2 - media*media) / n);      // Altrimenti, calcola e restituisce l'incertezza statistica
}

// Funzione per calcolare il chi quadrato
double chi2(unsigned int oss, int teo){
    return ((oss - teo)*(oss - teo))/static_cast<double>(teo);
}

// Funzione per campionare un angolo tra 0 e pi
double angle(Random& rnd){ 
    // Genera coordinate casuali x e y all'interno del semicerchio unitario superiore
    double x = rnd.Rannyu(-1, 1);
    double y = rnd.Rannyu();
    
    // Verifica se il punto è all'interno del cerchio unitario superiore
    if(x*x + y*y <= 1) return acos(x/sqrt(x*x + y*y)); // Se sì, calcola e restituisce l'angolo corrispondente
        
    else return angle(rnd);                            // Se no, richiama ricorsivamente la funzione fino a trovare un punto all'interno del semicerchio unitario superiore
}
        
