#include "Funzioni.hpp"

double IntegralMC::Mean(const BaseFunction *f, const double inf, const double sup, const unsigned int step) {
    // Imposta i membri della classe con i valori passati come argomento
    m_inf = inf;
    m_sup = sup;
    m_step = step;
    
	double x = 0.;
    double sum = 0.;

    // Ciclo per eseguire l'integrazione Monte Carlo
    for(unsigned int i = 0; i < m_step; i++) {
        x = m_gen.Rannyu(m_inf, m_sup); // Genera un numero casuale nell'intervallo [m_inf, m_sup]
        sum += f->Eval(x); // Valuta la funzione nel punto generato e aggiorna la sum 
    }
    
    // Calcola il valore medio dell'integranda moltiplicato per l'ampiezza dell'intervallo
    return (sum / m_step)*(m_sup - m_inf);
}   

// Importance Sampling
double IntegralMC::Importance(const BaseFunction *p, const BaseFunction *g, const double inf, const double sup,
	                          const double max, const unsigned int step){
	// Imposta i membri della classe con i valori passati come argomento
	m_inf = inf;
    m_sup = sup;
    m_step = step;

	double x = 0.;
    double sum = 0.;

    // Ciclo per eseguire l'integrazione Monte Carlo
    for(unsigned int i = 0; i < m_step; i++) {
        x = A_R(p, max); // Genera un numero casuale nell'intervallo [m_inf, m_sup] campionando p(x)
        sum += g->Eval(x); // Valuta la funzione g nel punto generato e aggiorna la sum 
    }
    
    // Calcola il valore medio della funzione g, che è una stima dell'integrale di partenza
    return (sum / m_step);

}

// Metodo di campionamento Accept-Reject
double IntegralMC::A_R(const BaseFunction *f, const double max) {
	// Genera due numeri casuali, il primo tra gli estremi dell'asse x, il secondo tra 0,1
	double x = m_gen.Rannyu(m_inf, m_sup);
	double r = m_gen.Rannyu();
	// Accetta x solo se r < f(x)/max
	return (r < f->Eval(x)/max) ? x : A_R(f, max);
}
	
// Funzione per calcolare l'incertezza statistica
double err(double mean, double mean2, unsigned int n){ 
    // Al primo blocco pone l'incertezza uguale a zero
    if(n == 0) return 0;
    else { // Altrimenti, calcola e restituisce l'incertezza statistica
        return sqrt((mean2 - mean*mean) / n);
    }
}

