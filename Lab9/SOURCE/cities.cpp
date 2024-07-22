#include "cities.hpp"  // Include il file di intestazione che dichiara la classe City

// Metodo per aggiungere una città alla mappa
void City::add_city(unsigned int label, double x, double y) {
    _label = label;  // Imposta l'etichetta della città
    _label_to_position[label] = vec({x, y});  // Aggiunge la città alla mappa con la posizione (x, y)
}

// Metodo per ottenere la posizione di una città data l'etichetta
vec City::get_position(unsigned int label) const {
    // Cerca l'etichetta nella mappa
    auto it = _label_to_position.find(label);
    
    // Se l'etichetta è trovata, restituisce la posizione associata
    if (it != _label_to_position.end()) {
        return it->second;
    } else {
        // Se l'etichetta non è trovata, stampa un messaggio di errore e restituisce un vettore vuoto
        cerr << "Label not found" << endl;
        return vec();
    }
}
