#ifndef __Funzioni_hpp__
#define __Funzioni_hpp__

#include <cmath>
#include <vector>
#include "random.h" // Includi la libreria per i numeri casuali

using namespace std;

// Classe astratta per definire una funzione generica
class BaseFunction {
public:
    virtual double Eval(double x) const = 0; // Metodo virtuale puro per valutare la funzione in un punto
};

// Classe che rappresenta la funzione da integrare
class Integrand: public BaseFunction {
public:
    Integrand() {;} // Costruttore vuoto
    ~Integrand() {;} // Distruttore vuoto
    virtual double Eval(double x) const override {
        return (M_PI / 2) * cos((M_PI * x) / 2); // Implementazione del metodo per valutare la funzione
    };
};

// Classe che rappresenta la distribuzione di probabilità p(x) per l'importance sampling
class dist_p: public BaseFunction {
public:
	dist_p() {;} 
	~dist_p() {;} 
	virtual double Eval(double x) const override {
		return 1.5*(1 - x*x);
	};
};

// Classe che rappresenta la funzione g(x) tale che p(x)g(x) = f(x) per l'importance sampling
class f_g: public BaseFunction {
public:
	f_g()  {;}
	~f_g() {;}
	virtual double Eval(double x) const override {
		return (M_PI*cos(0.5*M_PI*x))/(3*(1 - x*x));
	};
};	

// Classe per eseguire l'integrazione mediante il metodo Monte Carlo
class IntegralMC {
public:
    IntegralMC(Random& rnd) {m_gen = rnd;} // Costruttore che prende un oggetto Random come argomento
    ~IntegralMC() {;} // Distruttore vuoto
    
    // Metodo della media
    double Mean(const BaseFunction *f, const double inf, const double sup, const unsigned int step);
	// Importance Sampling
	double Importance(const BaseFunction *p, const BaseFunction *g, const double inf, const double sup, 
					  const double max, const unsigned int step);
	// Metodo di campionamento Accept-Reject
	double A_R(const BaseFunction *f, const double max);
	// Prende in ingresso una funzione, e il massimo della funzione

private:
    Random m_gen; // Oggetto Random per generare numeri casuali
    double m_inf; // Limite inferiore di integrazione
    double m_sup; // Limite superiore di integrazione
    unsigned int m_step; // Numero di passi per l'integrazione
};

// Funzione per calcolare l'incertezza statistica
double err(double mean, double mean2, unsigned int n);
// Prende in ingresso media, media quadratica e numero di blocchi usati (n = N-1)

// Funzione inline per calcolare la distanza quadratica
template <typename T> inline double dist2(T x, T y, T z) {
    return (double(x) * x + double(y) * y + double(z) * z);
};

#endif // Fine del file delle intestazioni