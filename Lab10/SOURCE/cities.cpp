#include "cities.hpp"

// Aggiunge una città con una data etichetta e posizione alla mappa
void City::add_city(unsigned int label, double x, double y) {
    _label = label;  // Imposta l'etichetta della città
    _label_to_position[label] = vec({x, y});  // Aggiunge la città alla mappa
}

// Metodo per ottenere la posizione di una città data l'etichetta
vec City::get_position(unsigned int label) const {
    auto it = _label_to_position.find(label);  // Cerca l'etichetta nella mappa
    if (it != _label_to_position.end()) {  // Se l'etichetta è trovata
        return it->second;  // Restituisce la posizione della città
    } else {
        cerr << "Label not found" << endl;  // Errore se l'etichetta non è trovata
        return vec();  // Restituisce un vettore vuoto
    }
}

// Serializza i dati delle città in un vettore di doubles
vec City::serialize() const {
    vec data;  // Vettore per memorizzare i dati serializzati
    for (const auto& entry : _label_to_position) {  // Itera su tutte le città nella mappa
        // Aggiunge l'etichetta della città al vettore
        data = join_vert(data, vec({static_cast<double>(entry.first)}));
        // Aggiunge la coordinata x della città al vettore
        data = join_vert(data, vec({entry.second(0)}));
        // Aggiunge la coordinata y della città al vettore
        data = join_vert(data, vec({entry.second(1)}));
    }
    return data;  // Restituisce il vettore contenente i dati serializzati
}

// Deserializza i dati delle città da un vettore di doubles
void City::deserialize(const vec& data) {
    _label_to_position.clear();  // Pulisce la mappa esistente
    for (unsigned int i = 0; i < data.size(); i += 3) {  // Itera sui dati in passi di 3 (etichetta, x, y)
        unsigned int label = static_cast<unsigned int>(data(i));  // Estrae l'etichetta
        double x = data(i + 1);  // Estrae la coordinata x
        double y = data(i + 2);  // Estrae la coordinata y
        // Aggiunge la città alla mappa
        _label_to_position[label] = vec({x, y});
    }
}
