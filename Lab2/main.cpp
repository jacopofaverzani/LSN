// Includi librerie e file header necessari
#include "Funzioni.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

int main(int argc, char *argv[]){

    Random rnd; // Oggetto della classe Random per generare numeri casuali
    rnd.initialize();

    unsigned int M = 10000; // Numero di punti
    unsigned int N = 100; // Numero di blocchi
    unsigned int L = M/N; // Numero di step per blocco

    // Stima dell'integrale I con metodo Monte Carlo campionando con distribuzione uniforme
    // Vettori per memorizzare le medie e le medie quadratiche progressive
    vector<double> mean(N);
    vector<double> mean2(N);
    
    // Estremi di integrazione
    double inf = 0.;
    double sup = 1.;
    
    IntegralMC integral(rnd);   // Oggetto per calcolare l'integrale Monte Carlo
    Integrand *f = new Integrand();  // Funzione da integrare

    // Loop sui blocchi
    for(unsigned int i = 0; i < N; i++){
        mean[i] = 0.;
        // Loop sugli step all'interno del blocco
        for(unsigned int j = 0; j < L; j++){
            // Calcola l'integral con 1000 step
            double value = integral.Mean(f, inf, sup, 1000);
            mean[i] += value;
        }
        mean[i] /= L; // Calcola la media per il blocco corrente
        mean2[i] = mean[i]*mean[i]; // Calcola la media quadratica per il blocco corrente
    }

    ofstream out_unif("unif.txt");
    if(out_unif.is_open()){ 
        out_unif << "mean" << setw(20) << "err" << endl;
        // Calcola le medie progressive e le incertezze statistiche
        for (unsigned int i = 0; i < N; i++){
            double sum_prog = 0.;
            double sum_prog2 = 0.;
            double err_prog = 0.;
        
            for(unsigned int j = 0; j < i+1; j++){
                sum_prog += mean[j];   // Somma delle medie progressive fino al blocco corrente
                sum_prog2 += mean2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
            }
            sum_prog /= (i+1);     // Media progressiva cumulativa
            sum_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
            err_prog = err(sum_prog, sum_prog2, i); // Incertezza statistica
        
            // Scrive su file le medie progressive e le relative incertezze
            out_unif << sum_prog << setw(20) << err_prog << endl;
        }
    }else{
        cerr << "Problema nell'apertura del file unif.txt" << endl;
    }
    out_unif.close(); // Chiude il file di output

    // Stima dell'integrale I con metodo Monte Carlo utilizzando importance sampling
    // Definisco la distribuzione p(x) tale che f(x) = p(x)g(x)
    dist_p *p = new dist_p();
    f_g *g = new f_g();

    double max = 1.5; // Massimo della distribuzione p(x) nell'intervallo [0,1]

    // Loop sui blocchi
    for(unsigned int i = 0; i < N; i++){
        mean[i] = 0.;
        // Loop sugli step all'interno del blocco
        for(unsigned int j = 0; j < L; j++){
            // Calcola l'integrale con 1000 step
            double value = integral.Importance(p, g, inf, sup, max, 1000);
            mean[i] += value;
        }
        mean[i] /= L; // Calcola la media il blocco corrente
        mean2[i] = mean[i]*mean[i]; // Calcola la media quadratica per il blocco corrente
    }
    
    ofstream out_imp("importance.txt");
    if(out_imp.is_open()){ 
        out_imp << "mean" << setw(20) << "err" << endl;
        // Calcola le medie progressive e le incertezze statistiche
        for (unsigned int i = 0; i < N; i++){
            double sum_prog = 0.;
            double sum_prog2 = 0.;
            double err_prog = 0.;
        
            for(unsigned int j = 0; j < i+1; j++){
                sum_prog += mean[j];   // Somma delle medie progressive fino al blocco corrente
                sum_prog2 += mean2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
            }
            sum_prog /= (i+1);     // Media progressiva cumulativa
            sum_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
            err_prog = err(sum_prog, sum_prog2, i); // Incertezza statistica
        
            // Scrive su file le medie progressive e le relative incertezze
            out_imp << sum_prog << setw(20) << err_prog << endl;
        }
    }else{
        cerr << "Problema nell'apertura del file importance.txt" << endl;
    }
    out_imp.close(); // Chiude il file di output

    // Stima della radice della distanza quadratica media per Random Walk 3D a distanza fissata
    unsigned int len = 100; // Lunghezza massima del Random Walk

    // Apre il file per salvare i risultati del Random Walk su reticolo
    ofstream rw_r("rw_ret.txt");
    if (rw_r.is_open()) {
        rw_r << "dist" << setw(20) << "err" << endl;
        // Ciclo sulla lunghezza del Random Walk
        for (unsigned int length = 1; length <= len; length++) {
            double dist = 0.;
            double dist_2 = 0.;

            // Ciclo sui blocchi
            for (unsigned int i = 0; i < N; i++) {
                double sum_dist2 = 0.; // Somma di un blocco

                // Ciclo sugli step all'interno di un blocco
                for (unsigned int j = 0; j < L; j++) {
                    vector<int> pos(3, 0); // Posizione iniziale del Random Walk

                    // Esegue il Random Walk
                    for (unsigned int l = 0; l < length; l++) {
                        unsigned int dir = static_cast<unsigned int>(rnd.Rannyu() * 3); // Sceglie direzione casuale
                        int sense = (rnd.Rannyu() < 0.5) ? -1 : 1; // Sceglie senso casuale
                        pos[dir] += sense;
                    }
                    // Distanza quadratica alla fine del blocco
                    double value_dist2 = dist2(pos[0], pos[1], pos[2]); 
                    sum_dist2 += value_dist2;
                }
                sum_dist2 = sqrt(sum_dist2 / L);
                dist += sum_dist2;
                dist_2 += sum_dist2 * sum_dist2;
            }
            dist /= len;
            dist_2 /= len;
            double error = err(dist, dist_2, N - 1);
            rw_r << dist << setw(20) << error << endl;
        }
    } else {
        cerr << "Problem to open rw_ret.txt" << endl;
    }
    rw_r.close(); 

    // Apre il file per salvare i risultati del Random Walk nel continuo
    ofstream rw_c("rw_cont.txt");
    if (rw_c.is_open()) {
        rw_c << "dist" << setw(20) << "err" << endl;
        // Ciclo sulla lunghezza del Random Walk
        for (unsigned int length = 1; length <= len; length++) {
            double dist = 0.;
            double dist_2 = 0.;

            // Ciclo sui blocchi
            for (unsigned int i = 0; i < N; i++) {
                double sum_dist2 = 0.; // Somma di un blocco

                // Ciclo sugli step all'interno di un blocco
                for (unsigned int j = 0; j < L; j++) {
                    vector<double> pos(3, 0.0); // Posizione iniziale del Random Walk

                    // Esegue il Random Walk
                    for (unsigned int l = 0; l < length; l++) {
                        // Scelgo una direzione casuale in coordinate sferiche
                        double theta = rnd.Rannyu(0, M_PI); 
                        double phi = rnd.Rannyu(0, 2 * M_PI);
                        pos[0] += sin(theta) * cos(phi);
                        pos[1] += sin(theta) * sin(phi);
                        pos[2] += cos(theta);
                    }
                    // Distanza quadratica alla fine del blocco
                    double value_dist2 = dist2(pos[0], pos[1], pos[2]); 
                    sum_dist2 += value_dist2;
                }   
                sum_dist2 = sqrt(sum_dist2 / L);
                dist += sum_dist2;
                dist_2 += sum_dist2 * sum_dist2;
            }
            dist /= len;
            dist_2 /= len;
            double errore = err(dist, dist_2, N - 1);
            rw_c << dist << setw(20) << errore << endl;
        }
    } else {
        cerr << "Problem to open rw_cont.txt" << endl;
    }
    rw_c.close();

    rnd.SaveSeed();
    return 0;
}
