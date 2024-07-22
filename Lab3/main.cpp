// Includi librerie e file header necessari
#include "Funzioni.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

using namespace std;

int main(int argc, char *argv[]){
    
    Random rnd; // Oggetto della classe Random per generare numeri casuali
    rnd.initialize();

    unsigned int M = 10000; // Numero di punti
    unsigned int N = 100; // Numero di blocchi
    unsigned int L = M/N; // Numero di step per blocco

    // Vettori per memorizzare le medie e le medie quadratiche progressive
    vector<double> mean_c(N); // Call
    vector<double> mean2_c(N);
    vector<double> mean_p(N); // Put
    vector<double> mean2_p(N);

    // Parametri
    double r = 0.1; // Fattore di rischio
    double s = 0.25; // Volatilità
    double S0 = 100.; // Prezzo dell'asset a t = 0
    double T = 1.; // Tempo di scadenza
    double K = 100.; // Prezzo di esercizio

    // Calcola al tempo t = 0 il prezzo di una opzione "call" e "put" Europea tramite metodo Monte Carlo 
    // campionando direttamente il prezzo finale dell'asset S(T) per un GBM(r, s^2) 
    // Loop sui blocchi
    for(unsigned int i = 0; i < N; i++){
        mean_c[i] = 0.;
        mean_p[i] = 0.;
        // Loop sugli step all'interno del blocco
        for(unsigned int j = 0; j < L; j++){
            // Campiona direttamente il prezzo finale S(T) per un GBM(r, s^2) 
            double S = European_dir(r, s, S0, T, rnd);
            mean_c[i] += exp(-r*T)*max(0., S - K);
            mean_p[i] += exp(-r*T)*max(0., K - S);
        }
        mean_c[i] /= L; // Calcola la media per il blocco corrente
        mean2_c[i] = mean_c[i]*mean_c[i]; // Calcola la media quadratica per il blocco corrente
        mean_p[i] /= L; 
        mean2_p[i] = mean_p[i]*mean_p[i];
    }

    ofstream out_dir("diretto.txt");
    if(out_dir.is_open()){ 
        out_dir << "mean_c" << setw(20) << "err_c" << setw(20) << "mean_p" << setw(20) << "err_p" << endl;
        // Calcola le medie progressive e le incertezze statistiche
        for (unsigned int i = 0; i < N; i++){
            double sum_prog_c = 0.;
            double sum_prog2_c = 0.;
            double err_prog_c = 0.;
            double sum_prog_p = 0.;
            double sum_prog2_p = 0.;
            double err_prog_p = 0.;
        
            for(unsigned int j = 0; j < i+1; j++){
                sum_prog_c += mean_c[j];   // Somma delle medie progressive fino al blocco corrente
                sum_prog2_c += mean2_c[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
                sum_prog_p += mean_p[j];   
                sum_prog2_p += mean2_p[j];
            }
            sum_prog_c /= (i+1);     // Media progressiva cumulativa
            sum_prog2_c /= (i+1);    // Media quadratica progressiva cumulativa
            err_prog_c = err(sum_prog_c, sum_prog2_c, i); // Incertezza statistica
            sum_prog_p /= (i+1);
            sum_prog2_p /= (i+1);    
            err_prog_p = err(sum_prog_p, sum_prog2_p, i); 
        
            // Scrive su file le medie progressive e le relative incertezze
            out_dir << sum_prog_c << setw(20) << err_prog_c << setw(20) << sum_prog_p << setw(20) << err_prog_p << endl;
        }
    }else{
        cerr << "Problema nell'apertura del file diretto.txt" << endl;
    }
    out_dir.close(); // Chiude il file di output

    // Calcola al tempo t = 0 il prezzo di una opzione "call" e "put" Europea tramite metodo Monte Carlo 
    // campionando il percorso GBM(r, s^2) discreto  dell prezzo dell'asset dividendo [0,T] in 100 intervalli di tempo 

    // Loop sui blocchi
    for(unsigned int i = 0; i < N; i++){
        mean_c[i] = 0.;
        mean2_p[i] = 0.;
        // Loop sugli step all'interno del blocco
        for(unsigned int j = 0; j < L; j++){
            // Campiona direttamente il prezzo finale S(T) per un GBM(r, s^2) 
            double S = European_step(r, s, S0, T, rnd);
            mean_c[i] += exp(-r*T)*max(0., S - K);
            mean_p[i] += exp(-r*T)*max(0., K - S);
        }
        mean_c[i] /= L; // Calcola la media il blocco corrente
        mean2_c[i] = mean_c[i]*mean_c[i]; // Calcola la media quadratica per il blocco corrente
        mean_p[i] /= L; 
        mean2_p[i] = mean_p[i]*mean_p[i];
    }

    ofstream out_step("percorso.txt");
    if(out_step.is_open()){
        out_step << "mean_c" << setw(20) << "err_c" << setw(20) << "mean_p" << setw(20) << "err_p" << endl;
        // Calcola le medie progressive e le incertezze statistiche
        for (unsigned int i = 0; i < N; i++){
            double sum_prog_c = 0.;
            double sum_prog2_c = 0.;
            double err_prog_c = 0.;
            double sum_prog_p = 0.;
            double sum_prog2_p = 0.;
            double err_prog_p = 0.;
        
            for(unsigned int j = 0; j < i+1; j++){
                sum_prog_c += mean_c[j];   // Somma delle medie progressive fino al blocco corrente
                sum_prog2_c += mean2_c[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
                sum_prog_p += mean_p[j];   
                sum_prog2_p += mean2_p[j];
            }
            sum_prog_c /= (i+1);     // Media progressiva cumulativa
            sum_prog2_c /= (i+1);    // Media quadratica progressiva cumulativa
            err_prog_c = err(sum_prog_c, sum_prog2_c, i); // Incertezza statistica
            sum_prog_p /= (i+1);
            sum_prog2_p /= (i+1);    
            err_prog_p = err(sum_prog_p, sum_prog2_p, i); 
        
            // Scrive su file le medie progressive e le relative incertezze
            out_step << sum_prog_c << setw(20) << err_prog_c << setw(20) << sum_prog_p << setw(20) << err_prog_p << endl;
        }
    }else{
        cerr << "Problema nell'apertura del file percorso.txt" << endl;
    }
    out_step.close(); // Chiude il file di output

    rnd.SaveSeed();
    return 0;
}