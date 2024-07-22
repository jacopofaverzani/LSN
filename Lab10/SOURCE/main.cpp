#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "genetic.hpp"

int main(int argc, char *argv[]) {
    // Inizializza MPI
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Ottiene il rango del processo corrente
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Ottiene il numero totale di processi

    Genetic GA;
    GA.initialize(rank, size);  // Inizializza l'algoritmo genetico con il rango e la dimensione della comunicazione MPI
    
    // Ottiene le informazioni dalla classe Genetic
    unsigned int num_cities = GA.num_cities();
    unsigned int population = GA.population_size();
    unsigned int half_pop = GA.half_population();
    unsigned int generations = GA.generations();
    unsigned int mig_freq = GA.mig_freq();

    double average_local = 0.0;
    double average_best;

    struct {
        double value;  // Valore della funzione di perdita
        int rank;      // Rango del processo che ha ottenuto il valore
    } local_loss, best_loss;

    local_loss.rank = rank;  // Imposta il rango locale

    // Calcola la funzione di perdita per tutti i percorsi iniziali
    for (unsigned int i = 0; i < population; i++) {
        GA.loss_function(i);
    }
    
    GA.order();  // Ordina i percorsi in base alla loro perdita
    local_loss.value = GA.get_loss(0);  // Salva la perdita del miglior percorso locale

    // Scrive la perdita minima globale nel file di output
    std::ofstream out_loss("../OUTPUT/loss.out", std::ios::app);
    MPI_Reduce(&local_loss, &best_loss, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        out_loss << "0" << std::setw(22) << best_loss.value << std::endl;
    }

    // Scrive la perdita media globale nel file di output
    std::ofstream out_ave("../OUTPUT/ave_loss.out", std::ios::app);
    for (unsigned int i = 0; i < half_pop; i++) {
        average_local += GA.get_loss(i);
    }
    average_local /= half_pop;
    MPI_Reduce(&average_local, &average_best, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        out_ave << "0" << std::setw(22) << average_best << std::endl;
    }

    // Esegue l'algoritmo genetico per un certo numero di generazioni
    for (unsigned int i = 1; i < generations + 1; i++) {
        GA.new_population();  // Crea una nuova generazione di percorsi
        for (unsigned int j = 0; j < population; j++) {
            GA.loss_function(j);  // Calcola la funzione di perdita per tutti i percorsi
        }
        GA.order();  // Ordina i percorsi in base alla loro perdita
        local_loss.value = GA.get_loss(0);  // Salva la perdita del miglior percorso locale
        MPI_Reduce(&local_loss, &best_loss, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            out_loss << i << std::setw(22) << best_loss.value << std::endl;
        }

        average_local = 0.0;
        for (unsigned int j = 0; j < half_pop; j++) {
            average_local += GA.get_loss(j);
        }
        average_local /= half_pop;
        MPI_Reduce(&average_local, &average_best, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            out_ave << i << std::setw(22) << average_best << std::endl;
        }

        // Esegue la migrazione se è il turno di migrare
        if (i % mig_freq == 0) {
            GA.migrate();
            GA.order();  // Ordina i percorsi dopo la migrazione
        }
    }

    out_loss.close();
    out_ave.close();

    // Broadcast del miglior percorso a tutti i processi
    MPI_Bcast(&best_loss.rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Scrive le posizioni delle città del miglior percorso globale nel file di output
    if (rank == best_loss.rank) {
        std::ofstream out_pos("../OUTPUT/pos.out", std::ios::app);
        for (unsigned int i = 0; i < num_cities + 1; i++) {
            unsigned int label = GA.get_city_label(0, i);
            out_pos << label << std::setw(22) << GA.get_city_coordinates(label)(0) << std::setw(22) << GA.get_city_coordinates(label)(1) << std::endl;
        }
        out_pos.close();
    }
    
    MPI_Finalize();  // Finalizza MPI
    return 0;
}
