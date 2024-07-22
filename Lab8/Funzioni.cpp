#include "Funzioni.hpp"

// Calcola l'energia E(x, u, s) in funzione di x, u, e s
double VMC::E(double x, double u, double s) const {
    return (-0.5 * (exp(-pow(x - u, 2) / (2 * s * s)) / (s * s) * (pow((x - u) / s, 2) - 1)
             + exp(-pow(x + u, 2) / (2 * s * s)) / (s * s) * (pow((x + u) / s, 2) - 1))) /
           (exp(-pow(x - u, 2) / (2 * s * s)) + exp(-pow(x + u, 2) / (2 * s * s))) + pow(x, 4) - 2.5 * x * x;
}

// Calcola la probabilità p(x, u, s) in funzione di x, u, e s
double VMC::p(double x, double u, double s) const {
    return pow(exp(-pow(x - u, 2) / (2 * s * s)) + exp(-pow(x + u, 2) / (2 * s * s)), 2);
}

// Metodo per eseguire l'algoritmo di Metropolis con distribuzione uniforme
double VMC::Metropolis_unif(unsigned int n_blocchi, unsigned int step, double dim_step, double u, double s, bool SA) {
    double x, x_prop;
    double somma;
    vector<double> media(n_blocchi);
    vector<double> media_2(n_blocchi);
    double media_prog, media_prog_2, err_prog;
    unsigned int L = int(step / n_blocchi);
    unsigned int passi_acc;

    // Se SA è true, esegue l'algoritmo di Simulated Annealing
    if (SA == true) {
        ofstream out_acc("AccettazioneMetropolis_SA.out");
        if (out_acc.is_open() == false) cerr << "Errore: non riesco ad aprire AccettazioneMetropolis_SA.out";

        for (unsigned int i = 0; i < n_blocchi; i++) {
            passi_acc = 0;
            somma = 0.0;
            x = 0.0;

            for (unsigned int j = 0; j < L; j++) {
                x_prop = x + gen.Rannyu(-dim_step, dim_step);
                if (gen.Rannyu() < min(1., p(x_prop, u, s) / p(x, u, s))) {
                    x = x_prop;
                    passi_acc += 1;
                }
                somma += E(x, u, s);
            }
            out_acc << static_cast<double>(passi_acc) / (L) << endl;
            media[i] = somma / L;
            media_2[i] = media[i] * media[i];
        }

        out_acc.close();

        for (unsigned int i = 0; i < n_blocchi; i++) {
            media_prog = 0.0;
            media_prog_2 = 0.0;

            for (unsigned int j = 0; j < i + 1; j++) {
                media_prog += media[j]; // Somma delle medie progressive fino al blocco corrente
                media_prog_2 += media_2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
            }
            media_prog /= (i + 1); // Media progressiva cumulativa
            media_prog_2 /= (i + 1); // Media quadratica progressiva cumulativa
            err_prog = errore(media_prog, media_prog_2, i); // Incertezza statistica
        }

        set_err(err_prog);
        return media_prog;
    }
    // Altrimenti, esegue l'algoritmo di Metropolis standard
    else {
        ofstream out_acc("AccettazioneMetropolis.out");
        if (out_acc.is_open() == false) cerr << "Errore: non riesco ad aprire AccettazioneMetropolis.out";

        ofstream out_ist("Istogramma.out");
        if (out_ist.is_open() == false) cerr << "Errore: non riesco ad aprire Istogramma.out";

        for (unsigned int i = 0; i < n_blocchi; i++) {
            passi_acc = 0;
            somma = 0.0;
            x = gen.Rannyu(-dim_step, dim_step);

            for (unsigned int j = 0; j < L; j++) {
                x_prop = x + gen.Rannyu(-dim_step, dim_step);
                if (gen.Rannyu() < min(1., p(x_prop, u, s) / p(x, u, s))) {
                    x = x_prop;
                    passi_acc += 1;
                }
                somma += E(x, u, s);
                out_ist << x << endl;
            }
            out_acc << static_cast<double>(passi_acc) / (L) << endl;
            media[i] = somma / L;
            media_2[i] = media[i] * media[i];
        }
        out_acc.close();
        out_ist.close();

        ofstream out_H("EnergiaMigliore.out");
        if (out_H.is_open() == false) cerr << "Errore: non riesco ad aprire EnergiaMigliore.out";
        out_H << "Blocco:" << setw(12) << "Energia:" << setw(12) << "Errore:" << endl;

        for (unsigned int i = 0; i < n_blocchi; i++) {
            media_prog = 0.0;
            media_prog_2 = 0.0;

            for (unsigned int j = 0; j < i + 1; j++) {
                media_prog += media[j]; // Somma delle medie progressive fino al blocco corrente
                media_prog_2 += media_2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
            }
            media_prog /= (i + 1); // Media progressiva cumulativa
            media_prog_2 /= (i + 1); // Media quadratica progressiva cumulativa
            err_prog = errore(media_prog, media_prog_2, i); // Incertezza statistica

            out_H << i << setw(20) << media_prog << setw(20) << err_prog << endl;
        }

        out_H.close();
        return 0;
    }
}

// Metodo per eseguire il Simulated Annealing
void VMC::Simulated_Annealing(unsigned int n_blocchi, unsigned int step, unsigned int dim_step,
                              double T_0, double T_fin, double scala,
                              vector<double>& Energie, vector<double>& Errori, vector<double>& Temperature,
                              vector<double>& u_list, vector<double>& s_list) {

    double u = get_u();
    double s = get_s();
    double E = Metropolis_unif(n_blocchi, step, dim_step, u, s, true);
    double err = get_err();
    double E_min = E;
    double u_min = u;
    double s_min = s;
    double T = T_0;
    double u_prop, s_prop, E_prop;

    // Ciclo di Simulated Annealing
    while (T > T_fin) {
        Temperature.push_back(T);
        u_list.push_back(u);
        s_list.push_back(s);
        u_prop = u + gen.Rannyu(-0.1, 0.1);
        s_prop = s + gen.Rannyu(-0.1, 0.1);
        if (s_prop <= 0) continue; // Evitare valori negativi o zero
        E_prop = Metropolis_unif(n_blocchi, step, dim_step, u_prop, s_prop, true);
        if (E_prop < E || exp(-(1 / T) * (E_prop - E)) < gen.Rannyu()) {
            E = E_prop;
            u = u_prop;
            s = s_prop;
            err = get_err();
        }
        if (E_prop < E_min) {
            E_min = E_prop;
            u_min = u_prop;
            s_min = s_prop;
        }
        T *= scala;
        Energie.push_back(E);
        Errori.push_back(err);
    }
    set_u(u_min);
    set_s(s_min);
    set_E(E_min);
}

// Funzione per calcolare l'incertezza statistica
double errore(double media, double media2, int n) {
    // Al primo blocco pone l'incertezza uguale a zero
    if (n == 0) {
        return 0;
    } else { // Altrimenti, calcola e restituisce l'incertezza statistica
        return sqrt((media2 - media * media) / n);
    }
}
