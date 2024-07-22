#include <vector>

#include "random.h"
#include "Funzioni.hpp"

int main(int argc, char *argv[]){
    
    Random rnd; // Oggetto della classe Random per generare numeri casuali
    rnd.initialize();
    
    unsigned int M = 1000000; // Numero di passi
    unsigned int N = 100; // Numero di blocchi
    unsigned int L = M/N; // Numero di passi per blocco

    // Vettori per memorizzare le medie e le medie quadratiche progressive
    vector<double> media(N);
    vector<double> media2(N);

    // Stato fondamentale e probabilità di transizione uniforme

    double posizione[3] = {0., 0., 0.}; // Posizione iniziale della particella
    unsigned int passi_accettati = 0; // Numero di passi accettati dall'algoritmo di Metropolis
    double dim_step = 1.2; // Lunghezza di ogni passo nell'algoritmo
    stato_fondamentale *s = new stato_fondamentale();

    for(unsigned int j = 0; j < L*2; j++){ // Tempo di equilibrio
        passi_accettati = Metropolis_unif(posizione, dim_step, rnd, s, passi_accettati);
    }
    passi_accettati = 0;

    ofstream out_p_f_u("pos_f_u.xyz");
    if(out_p_f_u.is_open() == false) cerr << "Errore: non riesco ad aprire pos_f_u.xyz";
    for(unsigned int i = 0; i < N; i++){ // Loop sui blocchi
        double r = 0.; // Distanza dall'origine
        for(unsigned int j = 0; j < L; j++){ // Loop sul numero di passi per blocco
            passi_accettati = Metropolis_unif(posizione, dim_step, rnd, s, passi_accettati);
            r += sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
            if(j%10 == 0){ // Ogni 10 passi salvo la posizione
                out_p_f_u << posizione[0] << "    " << posizione[1] << "    " << posizione[2] << endl;
            }
        }
        media[i] = r/L; // Calcola la media per il blocco corrente
        media2[i] = media[i]*media[i]; // Calcola la media quadratica per il blocco corrente
    }   
    
    out_p_f_u.close();

    ofstream out_r_f_u("r_f_u.txt");
    if(out_r_f_u.is_open() == false) cerr << "Errore: non riesco ad aprire r_f_u.txt";
    // Calcola le medie progressive e le incertezze statistiche
    for (unsigned int i = 0; i < N; i++){
        double somma_prog = 0.;
        double somma_prog2 = 0.;
        double err_prog = 0.;
        
        for(unsigned int j = 0; j < i + 1; j++){
            somma_prog += media[j];   // Somma delle medie progressive fino al blocco corrente
            somma_prog2 += media2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
        }
        somma_prog /= (i+1);     // Media progressiva cumulativa
        somma_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
        err_prog = errore(somma_prog, somma_prog2, i); // Incertezza statistica
        
        // Scrive su file le medie progressive e le relative incertezze
        out_r_f_u << somma_prog << ",    " << err_prog << endl;
    }

    out_r_f_u.close();

    double accettazione = static_cast<double>(passi_accettati)/M; // Calcolo il rapporto di accettazione dell'algoritmo
    ofstream out_acc("accettazione.txt");
    if(out_acc.is_open() == false) cerr << "Errore: non riesco ad aprire accettazione.txt";
    out_acc << "Accettazione per lo stato fondamentale, probabilità di transizione uniforme:    " << accettazione << endl;

    // Stato eccitato e probabilità di transizione uniforme

    for (unsigned int i = 0; i < 3; ++i) { 
        posizione[i] = 3.; // Nuova posizione di partenza
    }  
    dim_step = 3.; 
    stato_eccitato *p = new stato_eccitato();

    for(unsigned int j = 0; j < L; j++){ // Tempo di equilibrio
        passi_accettati = Metropolis_unif(posizione, dim_step, rnd, p, passi_accettati);
    }
    passi_accettati = 0;

    ofstream out_pos_e_u("pos_e_u.xyz");
    if(out_pos_e_u.is_open() == false) cerr << "Errore: non riesco ad aprire pos_e_u.xyz";
    for(unsigned int i = 0; i < N; i++){ // Loop sui blocchi
        double r = 0.;
        for(unsigned int j = 0; j < L; j++){ // Loop sul numero di passi per blocco
            passi_accettati = Metropolis_unif(posizione, dim_step, rnd, p, passi_accettati);
            r += sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
            if(j%10 == 0){
                out_pos_e_u << posizione[0] << "    " << posizione[1] << "    " << posizione[2] << endl;
            }
        }
        media[i] = r/L; // Calcola la media per il blocco corrente
        media2[i] = media[i]*media[i]; // Calcola la media quadratica per il blocco corrente
    }   
    
    out_pos_e_u.close();

    ofstream out_r_e_u("r_e_u.txt");
    if(out_r_e_u.is_open() == false) cerr << "Errore: non riesco ad aprire r_e_u.txt";
    // Calcola le medie progressive e le incertezze statistiche
    for (unsigned int i = 0; i < N; i++){
        double somma_prog = 0.;
        double somma_prog2 = 0.;
        double err_prog = 0.;
        
        for(unsigned int j = 0; j < i + 1; j++){
            somma_prog += media[j];   // Somma delle medie progressive fino al blocco corrente
            somma_prog2 += media2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
        }
        somma_prog /= (i+1);     // Media progressiva cumulativa
        somma_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
        err_prog = errore(somma_prog, somma_prog2, i); // Incertezza statistica
        
        // Scrive su file le medie progressive e le relative incertezze
        out_r_e_u << somma_prog << ",    " << err_prog << endl;
    }

    out_r_e_u.close(); // Chiude il file di output

    accettazione = static_cast<double>(passi_accettati)/M;
    out_acc << "Accettazione per lo stato eccitato, probabilità di transizione uniforme:    " << accettazione << endl;

    // Stato fondamentale e probabilità di transizione normale multivariata

    for (unsigned int i = 0; i < 3; ++i) { 
        posizione[i] = 0.; // Nuova posizione di partenza
    }  
    dim_step = 0.75; // Lunghezza di ogni passo nell'algoritmo

    for(unsigned int j = 0; j < L; j++){ // Tempo di equilibrio
        passi_accettati = Metropolis_multi(posizione, dim_step, rnd, s, passi_accettati);
    }
    passi_accettati = 0;

    ofstream out_p_f_m("pos_f_m.xyz");
    if(out_p_f_m.is_open() == false) cerr << "Errore: non riesco ad aprire pos_f_m.xyz";
    for(unsigned int i = 0; i < N; i++){ // Loop sui blocchi
        double r = 0.; // Distanza dall'origine
        for(unsigned int j = 0; j < L; j++){ // Loop sul numero di passi per blocco
            passi_accettati = Metropolis_multi(posizione, dim_step, rnd, s, passi_accettati);
            r += sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
            if(j%10 == 0){ // Ogni 10 passi salvo la posizione
                out_p_f_m << posizione[0] << "    " << posizione[1] << "    " << posizione[2] << endl;
            }
        }
        media[i] = r/L; // Calcola la media per il blocco corrente
        media2[i] = media[i]*media[i]; // Calcola la media quadratica per il blocco corrente
    }   
    
    out_p_f_m.close();

    ofstream out_r_f_m("r_f_m.txt");
    if(out_r_f_m.is_open() == false) cerr << "Errore: non riesco ad aprire r_f_m.txt";
    // Calcola le medie progressive e le incertezze statistiche
    for (unsigned int i = 0; i < N; i++){
        double somma_prog = 0.;
        double somma_prog2 = 0.;
        double err_prog = 0.;
        
        for(unsigned int j = 0; j < i + 1; j++){
            somma_prog += media[j];   // Somma delle medie progressive fino al blocco corrente
            somma_prog2 += media2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
        }
        somma_prog /= (i+1);     // Media progressiva cumulativa
        somma_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
        err_prog = errore(somma_prog, somma_prog2, i); // Incertezza statistica
        
        // Scrive su file le medie progressive e le relative incertezze
        out_r_f_m << somma_prog << ",    " << err_prog << endl;
    }

    out_r_f_m.close();

    accettazione = static_cast<double>(passi_accettati)/M; // Calcolo il rapporto di accettazione dell'algoritmo

    out_acc << "Accettazione per lo stato fondamentale, probabilità di transizione normale multivariata:    " << 
    accettazione << endl;

    // Stato eccitato e probabilità di transizione normale multivariata

    for (unsigned int i = 0; i < 3; ++i) { 
        posizione[i] = 3.; // Nuova posizione di partenza
    }  
    dim_step = 1.85; 

    for(unsigned int j = 0; j < L; j++){ // Tempo di equilibrio
        passi_accettati = Metropolis_multi(posizione, dim_step, rnd, p, passi_accettati);
    }
    passi_accettati = 0;

    ofstream out_pos_e_m("pos_e_m.xyz");
    if(out_pos_e_m.is_open() == false) cerr << "Errore: non riesco ad aprire pos_e_m.xyz";
    for(unsigned int i = 0; i < N; i++){ // Loop sui blocchi
        double r = 0.;
        for(unsigned int j = 0; j < L; j++){ // Loop sul numero di passi per blocco
            passi_accettati = Metropolis_multi(posizione, dim_step, rnd, p, passi_accettati);
            r += sqrt(posizione[0]*posizione[0] + posizione[1]*posizione[1] + posizione[2]*posizione[2]);
            if(j%10 == 0){
                out_pos_e_m << posizione[0] << "    " << posizione[1] << "    " << posizione[2] << endl;
            }
        }
        media[i] = r/L; // Calcola la media per il blocco corrente
        media2[i] = media[i]*media[i]; // Calcola la media quadratica per il blocco corrente
    }   
    
    out_pos_e_m.close();

    ofstream out_r_e_m("r_e_m.txt");
    if(out_r_e_m.is_open() == false) cerr << "Errore: non riesco ad aprire r_e_m.txt";
    // Calcola le medie progressive e le incertezze statistiche
    for (unsigned int i = 0; i < N; i++){
        double somma_prog = 0.;
        double somma_prog2 = 0.;
        double err_prog = 0.;
        
        for(unsigned int j = 0; j < i + 1; j++){
            somma_prog += media[j];   // Somma delle medie progressive fino al blocco corrente
            somma_prog2 += media2[j]; // Somma delle medie quadratiche progressive fino al blocco corrente
        }
        somma_prog /= (i+1);     // Media progressiva cumulativa
        somma_prog2 /= (i+1);    // Media quadratica progressiva cumulativa
        err_prog = errore(somma_prog, somma_prog2, i); // Incertezza statistica
        
        // Scrive su file le medie progressive e le relative incertezze
        out_r_e_m << somma_prog << ",    " << err_prog << endl;
    }

    out_r_e_m.close(); // Chiude il file di output

    accettazione = static_cast<double>(passi_accettati)/M;
    out_acc << "Accettazione per lo stato eccitato, probabilità di transizione normale multivariata:    " << accettazione << endl;

    out_acc.close();

    rnd.SaveSeed();
    
    return 0;
}