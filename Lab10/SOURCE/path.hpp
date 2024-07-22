#ifndef __path__
#define __path__

#include <armadillo>  // Include la libreria Armadillo per l'uso di vettori e matrici

using namespace std;  // Usa lo spazio dei nomi standard
using namespace arma;  // Usa lo spazio dei nomi Armadillo

// Classe che rappresenta un percorso nel problema del TSP
class Path {
private:
    vec _path;  // Vettore che memorizza il percorso delle città
    double _loss;  // Valore della funzione di costo del percorso

public:
    // Restituisce il vettore del percorso
    vec get_path() const { return _path; }

    // Imposta la dimensione del vettore del percorso
    void set_path_size(unsigned int size) { _path.set_size(size); }

    // Restituisce l'etichetta della città alla posizione specificata nel percorso
    unsigned int get_city_label(unsigned int pos) const { return _path(pos); }

    // Restituisce il valore della funzione di costo del percorso
    double get_loss() const { return _loss; }

    // Imposta l'etichetta della città alla posizione specificata nel percorso
    void set_city_label(unsigned int label, unsigned int pos) { _path(pos) = label; }

    // Imposta il valore della funzione di costo del percorso
    void set_loss(double loss) { _loss = loss; }
};

#endif  // Fine del blocco di inclusione condizionale
