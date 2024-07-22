#ifndef __genetic__    // Verifica se il simbolo __genetic__ non è già stato definito
#define __genetic__    // Definisce il simbolo __genetic__ per evitare inclusioni multiple

#include <fstream>         // Per la gestione dei file
#include <iomanip>        // Per la manipolazione dell'output
#include <cmath>          // Per funzioni matematiche
#include <algorithm>      // Per funzioni di algoritmo
#include <cassert>        // Per l'uso delle asserzioni
#include <unordered_set>  // Per l'uso degli insiemi non ordinati
#include "cities.hpp"     // Include la definizione della classe City
#include "path.hpp"       // Include la definizione della classe Path

// Classe che rappresenta l'algoritmo genetico
class Genetic {
private:
    unsigned int _sim_type;          // Tipo di simulazione
    unsigned int _num_cities;        // Numero di città
    unsigned int _population_size;   // Dimensione della popolazione
    unsigned int _generations;       // Numero di generazioni
    double _sel_exp;                 // Esponente di selezione
    vec _prob_mutations;             // Probabilità di mutazione
    double _prob_crossover;          // Probabilità di crossover
    City _city;                      // Oggetto di tipo City per gestire le città
    field<Path> _path;               // Campo per memorizzare i percorsi
    Random _rnd;                     // Oggetto per la generazione di numeri casuali

public:
    // Metodi getter per ottenere le variabili private
    unsigned int num_cities() const { return _num_cities; }
    unsigned int population_size() const { return _population_size; }
    unsigned int generations() const { return _generations; }
    vec get_city_coordinates(unsigned int label) const { return _city.get_position(label); }
    vec get_path(unsigned int i) const { return _path(i).get_path(); }
    unsigned int get_city_label(unsigned int path, unsigned int pos) const { return _path(path).get_city_label(pos); }
    double get_loss(unsigned int i) const { return _path(i).get_loss(); }

    // Dichiarazioni dei metodi pubblici della classe Genetic
    void initialize();          // Metodo per inizializzare l'algoritmo genetico
    void check_path(const vec path);  // Metodo per verificare la validità di un percorso
    void loss_function(unsigned int path); // Metodo per calcolare la funzione di perdita di un percorso
    void order();              // Metodo per ordinare la popolazione
    void selection(unsigned int& parent_1, unsigned int& parent_2); // Metodo per la selezione dei genitori
    void mutation_1(vec& path);  // Metodo di mutazione 1
    void mutation_2(vec& path);  // Metodo di mutazione 2
    void mutation_3(vec& path);  // Metodo di mutazione 3
    void mutation_4(vec& path);  // Metodo di mutazione 4
    void crossover(vec& parent_1, vec& parent_2); // Metodo di crossover
    void new_population();     // Metodo per creare una nuova popolazione
};

#endif // Fine del blocco di inclusione condizionale
