#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "random.h"

using namespace std;

// Classe che rappresenta la funzione d'onda della particella
class VMC {
private:
    double u, s, E_min, err;
    Random gen;

public:
    // Costruttore della classe VMC
    VMC(double media, double sigma, Random& rnd) {
        u = media;
        s = sigma;
        gen = rnd;
        E_min = 0.0;
        err = 0.0;
    }

    // Distruttore della classe VMC
    ~VMC() { ; }

    // Metodi getter per ottenere i valori delle variabili private
    double get_u() const { return u; }
    double get_s() const { return s; }
    double get_E() const { return E_min; }
    double get_err() const { return err; }

    // Metodi setter per impostare i valori delle variabili private
    void set_u(double new_u) { u = new_u; }
    void set_s(double new_s) { s = new_s; }
    void set_E(double new_E_min) { E_min = new_E_min; }
    void set_err(double new_err) { err = new_err; }

    // Metodo per calcolare l'energia
    double E(double x, double u, double s) const;

    // Metodo per calcolare la probabilità
    double p(double x, double u, double s) const;

    // Metodo per eseguire il Metropolis con distribuzione uniforme
    double Metropolis_unif(unsigned int n_blocchi, unsigned int step, double dim_step, double u, double s, bool SA);

    // Metodo per eseguire il Simulated Annealing
    void Simulated_Annealing(unsigned int n_blocchi, unsigned int step, unsigned int dim_step, 
                             double T_0, double T_fin, double scala,
                             vector<double>& Energie, vector<double>& Errori, vector<double>& Temperature,
                             vector<double>& u_list, vector<double>& s_list);
};

// Funzione per calcolare l'incertezza statistica
double errore(double media, double media2, int n) {
    // Prende in ingresso la media, la media quadratica e il numero di blocchi usati (n = N-1)
    return sqrt((media2 - media*media) / n);
}
