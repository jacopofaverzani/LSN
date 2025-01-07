#include "genetic.hpp"

// Inizializza l'oggetto Genetic in base al contenuto dei file di input nella directory ../INPUT/
void Genetic::initialize() {
    int p1, p2; // Numeri letti dal file ../INPUT/Primes per inizializzare il generatore di numeri casuali
    ifstream Primes("../INPUT/primes64001.in");
    Primes >> p1 >> p2;
    Primes.close();

    int seed[4]; // Seed per il generatore di numeri casuali
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    Seed.close();
    _rnd.SetRandom(seed, p1, p2); // Inizializza il generatore di numeri casuali

    ofstream output("../OUTPUT/output.out");
    if (!output.is_open()) {
        cerr << "Errore: non riesco ad aprire output.out" << endl;
    }

    ifstream input("../INPUT/input.dat"); // Inizio lettura del file di input
    if (!input.is_open()) {
        cerr << "Errore: non riesco ad aprire input.dat" << endl;
    }

    string property;

    // Lettura delle proprietà dal file di input
    while (!input.eof()) {
        input >> property;
        if (property == "TIPO_SIMULAZIONE") {
            input >> _sim_type;
            if (_sim_type == 0) {
                output << "Tipo simulazione: circonferenza" << endl;
            } else {
                output << "Tipo simulazione: quadrato" << endl;
            }
        }
        else if (property == "CITTA") {
            input >> _num_cities;
            output << "Numero di città: " << _num_cities << endl;
        }
        else if (property == "POPOLAZIONE") {
            input >> _population_size;
            _path.set_size(_population_size);
            output << "Dimensione della popolazione: " << _population_size << endl;
        }
        else if (property == "GENERAZIONI") {
            input >> _generations;
            output << "Numero di generazioni: " << _generations << endl;
        }
        else if (property == "PROBABILITA_MUTAZIONI") {
            _prob_mutations.set_size(4);
            input >> _prob_mutations(0) >> _prob_mutations(1) >> _prob_mutations(2) >> _prob_mutations(3);
            output << "Probabilità mutazione 1: " << _prob_mutations(0) << endl
                   << "Probabilità mutazione 2: " << _prob_mutations(1) << endl
                   << "Probabilità mutazione 3: " << _prob_mutations(2) << endl
                   << "Probabilità mutazione 4: " << _prob_mutations(3) << endl;
        }
        else if (property == "ESPONENTE_SELEZIONE") {
            input >> _sel_exp;
            output << "Esponente selezione: " << _sel_exp << endl;
        }
        else if (property == "PROBABILITA_CROSSOVER") {
            input >> _prob_crossover;
            output << "Probabilità crossover: " << _prob_crossover << endl;
        }
        else if (property == "FINE_INPUT") {
            break;
        } else {
            cerr << "Errore nella lettura dell'input." << endl;
        }
    }
    input.close();

    // Inizializza le città
    if (_sim_type == 0) { // Simulazione circonferenza
        for (unsigned int i = 0; i < _num_cities; i++) {
            double theta = _rnd.Rannyu(0, 2 * M_PI); // Angolo casuale
            _city.add_city(i, cos(theta), sin(theta)); // Aggiungi città con coordinate sulla circonferenza
        }
    } else { // Simulazione quadrato
        for (unsigned int i = 0; i < _num_cities; i++) {
            double x = _rnd.Rannyu(); // Coordinata x casuale
            double y = _rnd.Rannyu(); // Coordinata y casuale
            _city.add_city(i, x, y); // Aggiungi città con coordinate casuali
        }
    }

    // Inizializza la popolazione di percorsi
    for (unsigned int i = 0; i < _population_size; i++) {
        _path(i).set_path_size(_num_cities + 1);
        _path(i).set_city_label(0, 0); // Inizio del percorso
        _path(i).set_city_label(0, _num_cities); // Fine del percorso
        vec labels = arma::linspace<vec>(1, _num_cities - 1, _num_cities - 1); // Etichette delle città
        for (unsigned int j = 0; j < 10 * _num_cities; j++) {
            unsigned int index_1 = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1));
            unsigned int index_2 = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1));
            swap(labels(index_1), labels(index_2)); // Shuffle delle etichette
        }
        for (unsigned int j = 0; j < labels.size(); j++) {
            _path(i).set_city_label(static_cast<unsigned int>(labels(j)), j + 1);
        }
        check_path(_path(i).get_path()); // Verifica del percorso
    }

    // Creazione dei file di output per la loss e la media delle loss
    ofstream out_loss("../OUTPUT/loss.out");
    if (!out_loss.is_open()) {
        cerr << "Errore: non riesco ad aprire loss.out" << endl;
    }
    out_loss << "Generazione: " << setw(22) << "Migliore valore della funzione di costo: " << endl;
    out_loss.close();

    ofstream out_ave("../OUTPUT/ave_loss.out");
    if (!out_ave.is_open()) {
        cerr << "Errore: non riesco ad aprire ave_loss.out" << endl;
    }
    out_ave << "Generazione: " << setw(12) << "Media: " << endl;
    out_ave.close();

    ofstream out_pos("../OUTPUT/pos.out");
    if (!out_pos.is_open()) {
        cerr << "Errore: non riesco ad aprire pos.out" << endl;
    }
    out_pos << "Città: " << setw(12) << "x: " << setw(22) << "y: " << endl;
    out_pos.close();

    output << "Inizializzazione completata!" << endl;
    output.close();
}

// Verifica che il percorso sia valido
void Genetic::check_path(const vec path) {
    assert(path(0) == 0 && "Errore: la prima città del percorso non è quella iniziale.");
    assert(path(_num_cities) == 0 && "Errore: l'ultima città del percorso non è quella iniziale.");
    unordered_set<unsigned int> elements(path.memptr() + 1, path.memptr() + _num_cities);
    assert(elements.size() == _num_cities - 1 && "L'insieme deve contenere tutti i numeri da 1 a _num_cities - 1 esattamente una volta");
}

// Calcola la funzione di costo per un percorso
void Genetic::loss_function(unsigned int i) {
    double loss = 0.0;
    vec distance;
    distance.set_size(2);
    for (unsigned int j = 0; j < _num_cities; j++) {
        unsigned int label_1 = _path(i).get_city_label(j);
        unsigned int label_2 = _path(i).get_city_label(j + 1);
        distance(0) = _city.get_position(label_1)(0) - _city.get_position(label_2)(0);
        distance(1) = _city.get_position(label_1)(1) - _city.get_position(label_2)(1);
        loss += arma::dot(distance, distance); // Calcolo della distanza euclidea
    }
    _path(i).set_loss(loss); // Imposta la perdita per il percorso
}

// Ordina la popolazione di percorsi in base al valore della funzione di costo
void Genetic::order() {
    auto compare = [](const Path &a, const Path &b) {
        return a.get_loss() < b.get_loss(); // Ordina i percorsi in base al valore della funzione di costo
    };

    // Converti il field in un std::vector per usare std::sort
    vector<Path> path_vec(_path.size());
    for (unsigned int i = 0; i < _path.size(); i++) {
        path_vec[i] = _path(i);
    }

    // Ordina il vettore
    sort(path_vec.begin(), path_vec.end(), compare);

    // Copia il vettore ordinato di nuovo nel field
    for (unsigned int i = 0; i < _path.size(); i++) {
        _path(i) = path_vec[i];
    }
}

// Seleziona due genitori per la generazione successiva
void Genetic::selection(unsigned int& parent_1, unsigned int& parent_2) {
    parent_1 = static_cast<unsigned int>(_population_size * pow(_rnd.Rannyu(), _sel_exp));
    parent_2 = static_cast<unsigned int>(_population_size * pow(_rnd.Rannyu(), _sel_exp));
}

// Applicazione della mutazione 1: scambio di due città
void Genetic::mutation_1(vec& path) {
    if (_prob_mutations(0) > _rnd.Rannyu()) {
        unsigned int index_1 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        unsigned int index_2 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        // Assicurati che index_1 e index_2 siano diversi
        while (index_1 == index_2) {
            index_2 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        }
        unsigned int temp = path(index_1);
        path(index_1) = path(index_2);
        path(index_2) = temp;
        check_path(path); // Verifica del percorso
    }
}

// Applicazione della mutazione 2: inversione di un sotto-percorso
void Genetic::mutation_2(vec& path) {
    if (_prob_mutations(1) > _rnd.Rannyu()) {
        vec sub_path = path.subvec(1, _num_cities - 1);
        unsigned int n = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 2));
        unsigned int m = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        unsigned int start = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1 - m));

        for (int i = m - 1; i > -1; i--) {
            unsigned int temp = sub_path(start + i);
            sub_path(start + i) = sub_path((start + i + n) % (_num_cities - 1));
            sub_path((start + i + n) % (_num_cities - 1)) = temp;
        }
        for (unsigned int i = 1; i < _num_cities; i++) {
            path(i) = sub_path(i - 1);
        }
        check_path(path); // Verifica del percorso
    }
}

// Applicazione della mutazione 3: scambio di due sotto-percorsi
void Genetic::mutation_3(vec& path) {
    if (_prob_mutations(2) > _rnd.Rannyu()) {
        vec sub_path = path.subvec(1, _num_cities - 1);
        unsigned int m = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities / 2));
        unsigned int start_1 = static_cast<unsigned int>(_rnd.Rannyu(0, (_num_cities - 1 - m) / 2));
        unsigned int start_2 = static_cast<unsigned int>(_rnd.Rannyu((_num_cities - 1 - m) / 2, _num_cities - 1 - m));

        for (unsigned int i = 0; i < m; i++) {
            unsigned int temp = sub_path(start_1 + i);
            sub_path(start_1 + i) = sub_path(start_2 + i);
            sub_path(start_2 + i) = temp;
        }
        for (unsigned int i = 1; i < _num_cities; i++) {
            path(i) = sub_path(i - 1);
        }
        check_path(path); // Verifica del percorso
    }
}

// Applicazione della mutazione 4: inversione di un sotto-percorso
void Genetic::mutation_4(vec& path) {
    if (_prob_mutations(3) > _rnd.Rannyu()) {
        vec sub_path = path.subvec(1, _num_cities - 1);
        unsigned int m = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1));
        unsigned int start = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1 - m));

        for (unsigned int i = 0; i < static_cast<unsigned int>(m / 2); i++) {
            unsigned int temp = sub_path(start + i);
            sub_path(start + i) = sub_path(start + m - 1 - i);
            sub_path(start + m - 1 - i) = temp;
        }
        for (unsigned int i = 1; i < _num_cities; i++) {
            path(i) = sub_path(i - 1);
        }
        check_path(path); // Verifica del percorso
    }
}

// Applicazione del crossover tra due percorsi
void Genetic::crossover(vec& parent_1, vec& parent_2) {
    if (_prob_crossover > _rnd.Rannyu()) {
        unsigned int cut_index = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 2));

        vec son_1 = parent_1.subvec(1, cut_index);
        vec son_2 = parent_2.subvec(1, cut_index);

        unsigned int counter_1 = cut_index + 1;
        unsigned int counter_2 = cut_index + 1;

        // Riempie la parte rimanente di son_1 con città da parent_2
        for (unsigned int i = 1; i < _num_cities; i++) {
            if (!any(son_1 == parent_2(i))) {
                son_1.resize(counter_1);
                son_1(counter_1 - 1) = parent_2(i);
                counter_1++;
            }
        }

        // Riempie la parte rimanente di son_2 con città da parent_1
        for (unsigned int i = 1; i < _num_cities; i++) {
            if (!any(son_2 == parent_1(i))) {
                son_2.resize(counter_2);
                son_2(counter_2 - 1) = parent_1(i);
                counter_2++;
            }
        }

        // Aggiorna i genitori con i nuovi figli
        for (unsigned int i = 1; i < _num_cities; i++) {
            parent_1(i) = son_1(i - 1);
            parent_2(i) = son_2(i - 1);
        }

        // Verifica se i percorsi sono validi
        check_path(parent_1);
        check_path(parent_2);
    }
}

// Crea una nuova popolazione di percorsi
void Genetic::new_population() {
    field<Path> new_gen;
    new_gen.set_size(_population_size);
    for (unsigned int i = 0; i < _population_size / 2; i++) {
        unsigned int index_1, index_2;
        selection(index_1, index_2);
        vec parent_1 = _path(index_1).get_path();
        vec parent_2 = _path(index_2).get_path();
        crossover(parent_1, parent_2);

        mutation_1(parent_1);
        mutation_2(parent_1);
        mutation_3(parent_1);
        mutation_4(parent_1);

        mutation_1(parent_2);
        mutation_2(parent_2);
        mutation_3(parent_2);
        mutation_4(parent_2);

        new_gen(i).set_path_size(_num_cities + 1);
        new_gen(_population_size - 1 - i).set_path_size(_num_cities + 1);
        for (unsigned int j = 0; j < _num_cities + 1; j++) {
            new_gen(i).set_city_label(parent_1(j), j);
            new_gen(_population_size - 1 - i).set_city_label(parent_2(j), j);
        }
        check_path(new_gen(i).get_path());
        check_path(new_gen(_population_size - 1 - i).get_path());
    }
    // Aggiorna la popolazione con la nuova generazione
    for (unsigned int i = 0; i < _population_size; i++) {
        _path(i) = new_gen(i);
    }
}
