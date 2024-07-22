#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include "funzioni.hpp"

using namespace std;

int main (int argc, char *argv[]){
    
    Random rnd;                                           // Oggetto della classe Random per generare numeri casuali
    rnd.initialize();                                     // Inizializzo il generatore

    // Test del generatore utilizzando 10000 punti divisi in 100 blocchi
    unsigned int M = 10000;                               // Numero di punti
    unsigned int N = 100;                                 // Numero di blocchi
    unsigned int L = M/N;                                 // Numero di step per blocco
   
    // Vettori per memorizzare le medie e le medie quadratiche progressive per valor medio e deviazione standard
    vector<double> mean_m(N);
    vector<double> mean2_m(N);
    vector<double> mean_s(N);
    vector<double> mean2_s(N);
   
    for (unsigned int i = 0; i < N; i++){                 // Ciclo sui blocchi
        mean_m[i] = 0.;
        mean_s[i] = 0.;

        for (unsigned int j = 0; j < L; j++){             // Ciclo sui passi per blocco
            double x = rnd.Rannyu();
            mean_m[i] += x;
            mean_s[i] += (x - 0.5)*(x - 0.5);
        }
        mean_m[i] /= L;                                   // Media del valor medio del generatore al blocco i
        mean2_m[i] = mean_m[i]*mean_m[i];                 // Media quadratica del valor medio del generatore al blocco i
        mean_s[i] /= L;                                   // Media della varianza del generatore al blocco i
        mean2_s[i] = mean_s[i]*mean_s[i];                 // Media quadratica della varianza del generatore al blocco i
    }

    ofstream test("test.txt");
    if(test.is_open()){
        test << "mu" << setw(12) << "err_mu" << setw(12) << "sigma" << setw(12) << "err_sigma" << endl;
        for (unsigned int i = 0; i < N; i++){
            double sum_prog_m = 0.;
            double sum_prog2_m = 0.;
            double err_prog_m = 0.;
            double sum_prog_s = 0.;
            double sum_prog2_s = 0.;
            double err_prog_s = 0.;

            for (unsigned int j = 0; j < i + 1; j++){     // Sommo fino al blocco i compreso
                sum_prog_m += mean_m[j];
                sum_prog2_m += mean2_m[j];
                sum_prog_s += mean_s[j];
                sum_prog2_s += mean2_s[j];
            }
            sum_prog_m /= (i + 1);
            sum_prog2_m /= (i + 1);
            err_prog_m = err(sum_prog_m, sum_prog2_m, i); // Incertezza statistica progressiva al blocco i
            sum_prog_s /= (i + 1);
            sum_prog2_s /= (i + 1);
            err_prog_s = err(sum_prog_s, sum_prog2_s, i);

            test << sum_prog_m << setw(12) << setw(12) << err_prog_m << setw(12) << sum_prog_s << setw(12) << err_prog_s << endl;
        }
    } else cerr << "Problem to open test.txt" << endl;
    test.close();

    ofstream out_chi2("chi2.txt");
    if (out_chi2.is_open()) {
      out_chi2 << "chi2" << endl;  
      // Calcolo il chi quadro per 100 volte
      for (unsigned int i = 0; i < 100; i++) {
          double chi = 0.0;
          unsigned int cont;
          // Divido l'intervallo [0,1] in 100 blocchi
          for(unsigned int b = 0; b < N; b++){
            cont = 0;
            // Estraggo 10000 punti e conto quanti cadono in ciascun intervallo
            for(unsigned int j = 0; j < M; j++){
                double x = rnd.Rannyu();
                if (x >= b/static_cast<double>(N) && x < (b + 1)/static_cast<double>(N)) cont++;
            }
            chi += chi2(cont, int(M/N));
          }
          out_chi2 << chi << endl;
        } 
    } else cerr << "Problem to open chi2.txt" << endl;
    out_chi2.close();

    // Generazione di numeri casuali distribuiti uniformemente e secondo una distribuzione esponenziale e di Cauchy-Lorentz
    // per la verifica del Teorema del Limite Centrale
    // Scrittura dei risultati su file "unif.txt", "exp.txt" e "Lor.txt"
    // I risultati vengono generati per 10000 punti divisi in 100 blocchi e vengono calcolate le medie dei blocchi
    unsigned int ind[4] = {1, 2, 10, 100}; // Vettore di indici per il numero di variabili da sommare
    
    // Distribuzione uniforme
    ofstream out_unif("unif.txt");
    if(out_unif.is_open()){
        out_unif << "1" << setw(12) << "2" << setw(12) << "10" << setw(12) << "100" << endl;
        for(unsigned int i = 0; i < M; i++){
            for(unsigned int index : ind){
                double sum = 0.;
                for(unsigned int j = 0; j < index; j++){
                    double x = rnd.Rannyu();
                    sum += x;
                }
                sum /= index;
                if (index != 100) out_unif << sum << setw(12);
                else out_unif << sum << endl;
            }
        }
    } else cerr << "Problem to open unif.txt" << endl;
    out_unif.close();

    // Distribuzione esponenziale
    ofstream out_exp("exp.txt");
    if(out_exp.is_open()){
        out_exp << "1" << setw(12) << "2" << setw(12) << "10" << setw(12) << "100" << endl;
        for(unsigned int i = 0; i < M; i++){
            for(unsigned int index : ind){
                double sum = 0.;
                for(unsigned int j = 0; j < index; j++){
                    double x = rnd.exp(1.);
                    sum += x;
                }
                sum /= index;
                if (index != 100) out_exp << sum << setw(12);
                else out_exp << sum << endl;
                }
            }
        } else cerr << "Problem to open exp.txt" << endl;
    out_exp.close();

    // Distribuzione di Cauchy-Lorentz
    ofstream out_lor("Lor.txt");
    if(out_lor.is_open()){
        out_lor << "1" << setw(20) << "2" << setw(20) << "10" << setw(20) << "100" << endl;
        for(unsigned int i = 0; i < M; i++){
            for(unsigned int index : ind){
                double sum = 0.;
                for(unsigned int j = 0; j < index; j++){
                    double x = rnd.Lorentz(0., 1.);
                    sum += x;
                }
                sum /= index;
                if (index != 100) out_lor << sum << setw(20);
                else out_lor << sum << endl;
            }
        }
    } else cerr << "Problem to open Lor.txt" << endl;
    out_lor.close();

    // Simulazione dell'esperimento di Buffon per stimare il valore di π e scrittura dei risultati su file "Buffon.txt"
    // La lunghezza dell'ago è l, la distanza tra le linee è d
    vector<double> mean_pi(N);
    vector<double> mean2_pi(N);
    double l = 1.;
    double d = 2.;
    for(unsigned int i = 0; i < N; i++){
        mean_pi[i] = 0;
        for(unsigned int j = 0; j < L; j++){
            unsigned int N_int = 0; // Variabile per il numero di aghi che intersecano le linee
            double pi = 0;
            for(unsigned int n = 0; n < 1000; n++){
                double pos = rnd.Rannyu(0, d/2.); // Generazione della posizione del centro dell'ago
                double theta = angle(rnd);
                double len_eff = sin(theta)*l/2; // Rotazione dell'ago per un angolo casuale
                double x_prop = pos - len_eff;
                if(x_prop <= 0){ // Verifica dell'intersezione
                    N_int += 1;
                } 
            }
            pi = (2*l*1000)/(N_int*d);
            mean_pi[i] += pi;
        }
        mean_pi[i] /= L;
        mean2_pi[i] = mean_pi[i]*mean_pi[i];
    }

    // Scrittura su file delle stime progressive di π e delle relative incertezze statistiche
    ofstream out_buff("Buffon.txt");
    if(out_buff.is_open()){
        out_buff << "pi" << setw(12) << "err" << endl;
        for (unsigned int i = 0; i < N; i++){
            double sum_prog_pi = 0.;
            double sum_prog2_pi = 0.;
            double err_prog_pi = 0.;

            for(unsigned int j = 0; j < i + 1; j++){
                sum_prog_pi += mean_pi[j]; 
                sum_prog2_pi += mean2_pi[j];
            }
            sum_prog_pi /= (i+1);
            sum_prog2_pi /= (i+1);
            err_prog_pi = err(sum_prog_pi, sum_prog2_pi, i);

            out_buff << sum_prog_pi << setw(12) << err_prog_pi << endl;
        }
    } else cerr << "Problem to open Buffon.txt" << endl;
    out_buff.close();

    // Salva il seme finale del generatore di numeri casuali su "seed.out"
    rnd.SaveSeed();
    
    return 0;
}