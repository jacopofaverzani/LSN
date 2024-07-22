#include "Funzioni.hpp"
#include <algorithm>

// Funzione Metropolis_unif: implementa l'algoritmo Metropolis-Hastings con una distribuzione uniforme
unsigned int Metropolis_unif(double posizione[3], double dim_step, Random& rnd, const FunzioneBase_3D *p, 
                            unsigned int passi_accettati) {
    double pos_proposta[3]; // Posizione proposta per il prossimo passo
    double step[3]; // Passo di spostamento
    
    // Generazione della posizione proposta con distribuzione uniforme
    for(unsigned int i = 0; i < 3; i++) {
        step[i] = dim_step * rnd.Rannyu(-1, 1); 
        pos_proposta[i] = posizione[i] + step[i];
    }
    
    // Calcolo del fattore di accettazione
    double accettazione = min(1., p->Val(pos_proposta) / p->Val(posizione));
    
    // Decidere se accettare la nuova posizione
    if(rnd.Rannyu() < accettazione) {
        for(unsigned int i = 0; i < 3; i++) {
            posizione[i] = pos_proposta[i];
        }
        passi_accettati += 1; // Incremento del contatore dei passi accettati
    }
    return passi_accettati; // Restituisce il numero di passi accettati
}

// Funzione Metropolis_multi: implementa l'algoritmo Metropolis-Hastings con una distribuzione normale
unsigned int Metropolis_multi(double posizione[3], double dim_step, Random& rnd, const FunzioneBase_3D *p, 
                            unsigned int passi_accettati) {
    double pos_proposta[3]; // Posizione proposta per il prossimo passo
    double step[3]; // Passo di spostamento
    
    // Generazione della posizione proposta con distribuzione normale
    for(unsigned int i = 0; i < 3; i++) {
        step[i] = dim_step * rnd.Gauss(0., 1.); 
        pos_proposta[i] = posizione[i] + step[i];
    }
    
    // Calcolo del fattore di accettazione
    double accettazione = min(1., p->Val(pos_proposta) / p->Val(posizione));
    
    // Decidere se accettare la nuova posizione
    if(rnd.Rannyu() < accettazione) {
        for(unsigned int i = 0; i < 3; i++) {
            posizione[i] = pos_proposta[i];
        }
        passi_accettati += 1; // Incremento del contatore dei passi accettati
    }
    return passi_accettati; // Restituisce il numero di passi accettati
}                            

// Funzione errore: calcola l'incertezza statistica
double errore(double media, double media2, int n) { 
    // Al primo blocco pone l'incertezza uguale a zero
    if(n == 0) {
        return 0;
    } else { 
        // Altrimenti, calcola e restituisce l'incertezza statistica
        return sqrt((media2 - media * media) / n);
    }
}
