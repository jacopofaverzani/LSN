#include <iomanip>
#include "Funzioni.hpp"

int main(int argc, char *argv[]) {

    // Creazione di un oggetto Random per la generazione di numeri casuali
    Random rnd;
    rnd.initialize(); // Inizializzazione del generatore di numeri casuali
    
    unsigned int M = 20000; // Numero totale di passi per il Metropolis
    unsigned int N = 100; // Numero di blocchi per l'analisi

    double dim_step = 2.0; // Lunghezza di ciascun passo nell'algoritmo di Metropolis
    double u = 1.0, s = 1.0; // Valori iniziali per i parametri della funzione d'onda
    vector<double> Energie; // Vettore per memorizzare le energie ottenute
    vector<double> Errori; // Vettore per memorizzare gli errori associati alle energie
    vector<double> Temperature; // Vettore per memorizzare le temperature utilizzate
    vector<double> u_list; // Vettore per memorizzare i valori di u durante il Simulated Annealing
    vector<double> s_list; // Vettore per memorizzare i valori di s durante il Simulated Annealing

    // Creazione di un oggetto VMC con i parametri iniziali
    VMC *v = new VMC(u, s, rnd);

    double T_0 = 100.0; // Temperatura iniziale per il Simulated Annealing
    double T_fin = 1.0; // Temperatura finale per il Simulated Annealing
    double scala = 0.99; // Fattore di scalamento per la temperatura

    // Esecuzione del Simulated Annealing
    v->Simulated_Annealing(N, M, dim_step, T_0, T_fin, scala, Energie, Errori, Temperature, u_list, s_list);

    // Salvataggio dei risultati dell'annealing in un file di output
    ofstream out_H("H.out");
    if (out_H.is_open() == false) cerr << "Errore: non riesco ad aprire H.out";
    out_H << "Temperatura:" << setw(12) << "Energia:" << setw(16) << "Errore:" << setw(20) << "u:" << setw(20) << "s:" << endl;

    for (unsigned int i = 0; i < Energie.size(); i++) { // Loop sui blocchi
        out_H << Temperature[i] << setw(20) << Energie[i] << setw(20) << Errori[i] << setw(20) << u_list[i] << setw(20)
              << s_list[i] << endl;
    }
    out_H.close(); // Chiusura del file di output

    // Salvataggio dei migliori valori trovati in un file di output
    ofstream out_best("best.out");
    if (out_best.is_open() == false) cerr << "Errore: non riesco ad aprire best.out";
    out_best << "u best:" << setw(12) << "s best:" << setw(12) << "E_min:" << endl
    << v->get_u() << setw(12) << v->get_s() << setw(12) << v->get_E() << endl;
    out_best.close(); // Chiusura del file di output

    M = 100000; // Aggiornamento del numero di passi per il Metropolis
    // Esecuzione del Metropolis senza Simulated Annealing
    v->Metropolis_unif(N, M, dim_step, v->get_u(), v->get_s(), false);

    rnd.SaveSeed(); // Salvataggio del seme del generatore di numeri casuali
    return 0; // Fine del programma
}
