#ifndef __path__    // Verifica se il simbolo __path__ non è già stato definito
#define __path__    // Definisce il simbolo __path__ per evitare inclusioni multiple

#include <armadillo>    // Include la libreria Armadillo per la gestione di vettori e matrici

using namespace std;    // Usa lo spazio dei nomi standard per evitare di dover scrivere std:: ogni volta
using namespace arma;   // Usa lo spazio dei nomi di Armadillo

// Classe che rappresenta un percorso
class Path {
private:
    vec _path;    // Vettore che memorizza i percorsi delle città
    double _loss; // Valore che rappresenta il costo (o la perdita) del percorso

public:
    // Metodo per ottenere il vettore del percorso
    vec get_path() const { return _path; }

    // Metodo per impostare la dimensione del vettore del percorso
    void set_path_size(unsigned int size) { _path.set_size(size); }

    // Metodo per ottenere l'etichetta della città in una posizione specifica del percorso
    unsigned int get_city_label(unsigned int pos) const { return _path(pos); }

    // Metodo per ottenere il valore della loss del percorso
    double get_loss() const { return _loss; }

    // Metodo per impostare l'etichetta della città in una posizione specifica del percorso
    void set_city_label(unsigned int label, unsigned int pos) { _path(pos) = label; }

    // Metodo per impostare il valore della loss del percorso
    void set_loss(double loss) { _loss = loss; }
};

#endif // Fine del blocco di inclusione condizionale
