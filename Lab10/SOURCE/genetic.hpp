#ifndef __genetic__
#define __genetic__

#include <fstream>       // Per le operazioni di input/output su file
#include <iomanip>       // Per la formattazione dell'output
#include <cmath>         // Per funzioni matematiche avanzate
#include <algorithm>     // Per operazioni su container come sorting
#include <cassert>       // Per le asserzioni
#include <unordered_set> // Per l'uso di set non ordinati
#include <mpi.h>         // Per l'uso di MPI (Message Passing Interface)
#include "cities.hpp"    // Inclusione della classe City
#include "path.hpp"      // Inclusione della classe Path

using namespace std; // Usa lo spazio dei nomi standard
using namespace arma; // Usa lo spazio dei nomi Armadillo

// Classe che implementa un algoritmo genetico per il problema del TSP
class Genetic {
private:
    unsigned int _sim_type;          // Tipo di simulazione (es. circonferenza o quadrato)
    unsigned int _num_cities;        // Numero totale di città
    unsigned int _population_size;   // Dimensione della popolazione
    unsigned int _half_population;   // Dimensione della metà della popolazione
    unsigned int _sub_population;    // Dimensione della sottopopolazione per la migrazione
    unsigned int _generations;       // Numero di generazioni da eseguire
    unsigned int _mig_freq;          // Frequenza di migrazione tra le sottopopolazioni
    double _sel_exp;                 // Esponente per la selezione
    vec _prob_mutations;             // Probabilità di mutazione per diversi tipi di mutazione
    double _prob_crossover;          // Probabilità di crossover
    City _city;                      // Oggetto che rappresenta le città
    field<Path> _path;               // Campo di percorsi
    Random _rnd;                     // Generatore di numeri casuali
    unsigned int _rank;              // Rank del processo MPI
    unsigned int _size;              // Numero totale di processi MPI

public:
    // Funzioni di accesso per i membri privati
    unsigned int num_cities() const { return _num_cities; }
    unsigned int population_size() const { return _population_size; }
    unsigned int half_population() const { return _half_population; }
    unsigned int sub_population() const { return _sub_population; }
    unsigned int generations() const { return _generations; }
    unsigned int mig_freq() const { return _mig_freq; }
    vec get_city_coordinates(unsigned int label) const { return _city.get_position(label); }
    vec get_path(unsigned int i) const { return _path(i).get_path(); }
    unsigned int get_city_label(unsigned int path, unsigned int pos) const { return _path(path).get_city_label(pos); }
    double get_loss(unsigned int i) const { return _path(i).get_loss(); }

    // Funzioni per l'inizializzazione e la gestione dell'algoritmo genetico
    void initialize(int rank, int size);  // Inizializza il sistema con i parametri MPI
    void check_path(const vec path);      // Verifica la validità di un percorso
    void loss_function(unsigned int path); // Calcola la funzione di costo per un percorso
    void order();                        // Ordina i percorsi in base alla funzione di costo
    void selection(unsigned int& parent_1, unsigned int& parent_2); // Seleziona due genitori per la crossover
    void mutation_1(vec& path);           // Applicazione del primo tipo di mutazione
    void mutation_2(vec& path);           // Applicazione del secondo tipo di mutazione
    void mutation_3(vec& path);           // Applicazione del terzo tipo di mutazione
    void mutation_4(vec& path);           // Applicazione del quarto tipo di mutazione
    void crossover(vec& parent_1, vec& parent_2); // Applicazione della crossover tra due genitori
    void new_population();               // Crea una nuova popolazione a partire dalla popolazione attuale
    void migrate();                     // Gestisce la migrazione tra sottopopolazioni
};

#endif  // Fine del blocco di inclusione condizionale
