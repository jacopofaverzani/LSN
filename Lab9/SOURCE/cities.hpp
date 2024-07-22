#ifndef __cities__  // Verifica se il file non è già stato incluso
#define __cities__

#include <armadillo>   // Inclusione della libreria Armadillo per l'algebra lineare
#include <map>         // Inclusione della libreria standard per l'uso delle mappe
#include "random.hpp" // Inclusione di un file di intestazione per la gestione dei numeri casuali

using namespace std;  // Utilizza lo spazio dei nomi standard per evitare di dover scrivere std:: ogni volta
using namespace arma; // Utilizza lo spazio dei nomi Armadillo per evitare di dover scrivere arma:: ogni volta

// Classe che rappresenta una città con una posizione
class City {

private:
    // Mappa che associa un'etichetta (ID) a una posizione (vettore di coordinate)
    map<unsigned int, vec> _label_to_position;
    unsigned int _label; // Etichetta della città

public:
    // Aggiunge una città alla mappa con una data etichetta e posizione
    void add_city(unsigned int label, double x, double y);

    // Restituisce la posizione della città data un'etichetta
    vec get_position(unsigned int label) const;

    // Restituisce l'etichetta della città
    unsigned int get_label() const { return _label; }
};

#endif // Chiusura della direttiva di inclusione guard

