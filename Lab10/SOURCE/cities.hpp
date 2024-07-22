#ifndef __cities__
#define __cities__

#include <armadillo>
#include <map>
#include "random.hpp"

using namespace std;
using namespace arma;

class City {

private:
    // Mappa che associa un'etichetta (unsigned int) alla posizione della città (vector di double)
    map<unsigned int, vec> _label_to_position;
    
    // Etichetta della città (unsigned int)
    unsigned int _label;

public:
    // Aggiunge una città alla mappa con una data etichetta e posizione
    void add_city(unsigned int label, double x, double y);
    
    // Restituisce la posizione della città associata all'etichetta
    vec get_position(unsigned int label) const;

    // Restituisce l'etichetta della città
    unsigned int get_label() const { return _label; }
    
    // Serializza i dati delle città in un vettore di double
    vec serialize() const;
    
    // Deserializza i dati delle città da un vettore di double
    void deserialize(const vec& data);
};

#endif
