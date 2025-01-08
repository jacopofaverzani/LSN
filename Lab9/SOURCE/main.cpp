#include <iostream>
#include "genetic.hpp"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[]) {
    Genetic GA; // Crea un'istanza dell'oggetto Genetic
    GA.initialize(); // Inizializza l'oggetto Genetic con i parametri e le città

    // Ottieni i parametri della simulazione
    unsigned int num_cities = GA.num_cities(); // Numero di città
    unsigned int population = GA.population_size(); // Dimensione della popolazione
    unsigned int half_pop = population / 2; // Metà della popolazione
    unsigned int generations = GA.generations(); // Numero di generazioni
    double average = 0.0; // Media delle loss

    // Calcola la funzione di costo per ogni percorso iniziale
    for (unsigned int i = 0; i < population; i++) {
        GA.loss_function(i); // Calcola la loss per il percorso i
    }
    GA.order(); // Ordina la popolazione in base alla funzione di costo

    // Scrivi il valore iniziale della funzione di costo nel file
    ofstream out_loss("../OUTPUT/loss.out", ios::app);
    out_loss << "0" << setw(22) << GA.get_loss(0) << endl;

    // Calcola la media delle loss iniziali e scrivi nel file
    ofstream out_ave("../OUTPUT/ave_loss.out", ios::app);
    for (unsigned int i = 0; i < half_pop; i++) {
        average += GA.get_loss(i); // Somma le loss dei percorsi migliori
    }
    average /= half_pop; // Calcola la media
    out_ave << "0" << setw(22) << average << endl;

    // Ciclo per il numero di generazioni
    for (unsigned int i = 1; i < generations + 1; i++) {
        GA.new_population(); // Crea una nuova popolazione
        for (unsigned int j = 0; j < population; j++) {
            GA.loss_function(j); // Calcola la loss per ogni percorso
        }
        GA.order(); // Ordina la popolazione in base alla funzione di costo
        out_loss << i << setw(22) << GA.get_loss(0) << endl; // Scrivi il miglior valore della funzione di costo

        // Calcola la media delle loss e scrivi nel file
        average = 0.0;
        for (unsigned int j = 0; j < half_pop; j++) {
            average += GA.get_loss(j); // Somma le loss dei percorsi migliori
        }
        average /= half_pop; // Calcola la media
        out_ave << i << setw(22) << average << endl; // Scrivi la media delle loss
    }
    out_loss.close(); // Chiudi il file delle loss
    out_ave.close(); // Chiudi il file della media delle loss

    // Scrivi le coordinate delle città nel file di output
    ofstream out_pos("../OUTPUT/pos.out", ios::app);
    for (unsigned int i = 0; i < num_cities + 1; i++) {
        unsigned int label = GA.get_city_label(0, i); // Ottieni l'etichetta della città i
        out_pos << label << setw(22) << GA.get_city_coordinates(label)(0) << setw(22) << GA.get_city_coordinates(label)(1) << endl; // Scrivi l'etichetta e le coordinate
    }
    out_pos.close(); // Chiudi il file delle coordinate

    return 0; // Termina il programma
}
