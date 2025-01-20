#include "genetic.hpp"

void Genetic::initialize(int rank, int size) { 
    // Inizializza i membri della classe con i parametri di MPI
    _rank = rank;
    _size = size;

    vec p1_list, p2_list; // Vettori per i numeri primi letti dal file
    p1_list.set_size(_size);
    p2_list.set_size(_size);

    // Legge i numeri primi da ../INPUT/primes64001.in
    ifstream Primes("../INPUT/primes64001.in");
    if (!Primes.is_open()) {
        cerr << "Errore: non riesco ad aprire primes64001.in" << endl;
        return;
    }
    for (unsigned int i = 0; i < _size; i++) {
        Primes >> p1_list(i) >> p2_list(i);
    }
    Primes.close();

    // Legge il seme per il generatore di numeri casuali
    int seed[4];
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    Seed.close();
    _rnd.SetRandom(seed, p1_list(_rank), p2_list(_rank));

    // Crea e apre il file di output
    ofstream output("../OUTPUT/output.out");
    if (!output.is_open()) {
        cerr << "Errore: non riesco ad aprire output.out" << endl;
        return;
    }

    // Legge e processa il file di configurazione ../INPUT/input.dat
    ifstream input("../INPUT/input.dat");
    if (!input.is_open()) {
        cerr << "Errore: non riesco ad aprire input.dat" << endl;
        return;
    }
    string property;

    while (!input.eof()) {
        input >> property;
        if (property == "TIPO_SIMULAZIONE") {
            input >> _sim_type;
            if (_sim_type == 0) {
                output << "Tipo simulazione: circonferenza" << endl;
            } else if (_sim_type == 1) {
                output << "Tipo simulazione: quadrato" << endl;
            } else {
                output << "Tipo di simulazione: province" << endl;
            }
            output << "Numero di nodi: " << _size << endl;
        } 
        else if (property == "CITTA") {
            input >> _num_cities;
            output << "Numero di città: " << _num_cities << endl;
        } 
        else if (property == "POPOLAZIONE") {
            input >> _population_size;
            _half_population = _population_size / 2;
            _path.set_size(_population_size);
            output << "Dimensione della popolazione per ogni nodo: " << _population_size << endl;
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
        else if (property == "FREQUENZA_MIGRAZIONE") {
            input >> _mig_freq;
            output << "Frequenza migrazione: " << _mig_freq << endl;
        } 
        else if (property == "FINE_INPUT") {
            break;
        } 
        else {
            cerr << "Errore nella lettura dell'input." << endl;
        }
    }
    input.close();

    // Gestisce la simulazione in base al tipo specificato
    if (_sim_type == 0) { // Simulazione circonferenza
        if (rank == 0) {
            for (unsigned int i = 0; i < _num_cities; i++) {
                double theta = _rnd.Rannyu(0, 2 * M_PI);
                _city.add_city(i, cos(theta), sin(theta));
            }
        }
        // Serializza i dati delle città nel processo 0
        vec city_data;
        if (rank == 0) {
            city_data = _city.serialize();
        }

        // Trasmetti la dimensione dei dati a tutti i processi
        unsigned int data_size = city_data.size();
        MPI_Bcast(&data_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Alloca il buffer per i dati ricevuti negli altri processi
        if (rank != 0) {
            city_data.resize(data_size);
        }

        // Trasmetti i dati serializzati a tutti i processi
        MPI_Bcast(city_data.memptr(), data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Deserializza i dati nei processi diversi da 0
        if (rank != 0) {
            _city.deserialize(city_data);
        }
    } 
    else if (_sim_type == 1) { // Simulazione quadrato
        if (rank == 0) {
            for (unsigned int i = 0; i < _num_cities; i++) {
                double x = _rnd.Rannyu();
                double y = _rnd.Rannyu();
                _city.add_city(i, x, y);
            }
        }
        // Serializza i dati delle città nel processo 0
        vec city_data;
        if (rank == 0) {
            city_data = _city.serialize();
        }

        // Trasmetti la dimensione dei dati a tutti i processi
        unsigned int data_size = city_data.size();
        MPI_Bcast(&data_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Alloca il buffer per i dati ricevuti negli altri processi
        if (rank != 0) {
            city_data.resize(data_size);
        }

        // Trasmetti i dati serializzati a tutti i processi
        MPI_Bcast(city_data.memptr(), data_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Deserializza i dati nei processi diversi da 0
        if (rank != 0) {
            _city.deserialize(city_data);
        }
    } 
    else { // Simulazione province
        ifstream pos("../INPUT/cap_prov_ita.dat");
        if (!pos.is_open()) {
            cerr << "Errore: non riesco ad aprire cap_prov_ita.dat" << endl;
            return;
        }
        double x, y;
        for (unsigned int i = 0; i < _num_cities; i++) {
            pos >> x >> y;
            _city.add_city(i, x, y);
        }
        pos.close();
    }

    // Inizializza la popolazione
    for (unsigned int i = 0; i < _population_size; i++) {
        _path(i).set_path_size(_num_cities + 1);
        _path(i).set_city_label(0, 0);
        _path(i).set_city_label(0, _num_cities);
        vec labels = arma::linspace<vec>(1, _num_cities - 1, _num_cities - 1);
        for (unsigned int j = 0; j < 10 * _num_cities; j++) {
            unsigned int index_1 = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1));
            unsigned int index_2 = static_cast<unsigned int>(_rnd.Rannyu(0, _num_cities - 1));
            swap(labels(index_1), labels(index_2));
        }
        for (unsigned int j = 0; j < labels.size(); j++) {
            _path(i).set_city_label(static_cast<unsigned int>(labels(j)), j + 1);
        }
        check_path(_path(i).get_path());
    }

    // Crea e apre i file di output per le statistiche
    ofstream out_loss("../OUTPUT/loss.out");
    if (!out_loss.is_open()) {
        cerr << "Errore: non riesco ad aprire loss.out" << endl;
        return;
    }
    out_loss << "Generazione: " << setw(22) << "Migliore valore della funzione di costo: " << endl;
    out_loss.close();

    ofstream out_ave("../OUTPUT/ave_loss.out");
    if (!out_ave.is_open()) {
        cerr << "Errore: non riesco ad aprire ave_loss.out" << endl;
        return;
    }
    out_ave << "Generazione: " << setw(12) << "Media: " << endl;
    out_ave.close();

    ofstream out_pos("../OUTPUT/pos.out");
    if (!out_pos.is_open()) {
        cerr << "Errore: non riesco ad aprire pos.out" << endl;
        return;
    }
    out_pos << "Città: " << setw(12) << "x: " << setw(22) << "y: " << endl;
    out_pos.close();

    // Segnala la fine dell'inizializzazione
    output << "Inizializzazione completata!" << endl;
    output.close();
}


void Genetic::check_path(const vec path) {
    // Verifica che la prima e l'ultima città del percorso siano quelle iniziali
    assert(path(0) == 0 && "Errore: la prima città del percorso non è quella iniziale.");
    assert(path(_num_cities) == 0 && "Errore: l'ultima città del percorso non è quella iniziale.");

    // Verifica che tutte le città siano presenti una sola volta nel percorso
    unordered_set<unsigned int> elements(path.memptr() + 1, path.memptr() + _num_cities);
    assert(elements.size() == _num_cities - 1 && "L'insieme deve contenere tutti i numeri esattamente una volta");
}

void Genetic::loss_function(unsigned int i) {
    double loss = 0.0;
    vec distance;
    distance.set_size(2);

    // Calcola la funzione di perdita del percorso i-esimo
    for (unsigned int j = 0; j < _num_cities; j++) {
        unsigned int label_1 = _path(i).get_city_label(j);
        unsigned int label_2 = _path(i).get_city_label(j + 1);
        distance(0) = _city.get_position(label_1)(0) - _city.get_position(label_2)(0);
        distance(1) = _city.get_position(label_1)(1) - _city.get_position(label_2)(1);
        loss += arma::dot(distance, distance); // Quadrato della distanza euclidea
    }

    _path(i).set_loss(loss); // Imposta la perdita calcolata per il percorso i-esimo
}

void Genetic::order() {
    auto compare = [](const Path &a, const Path &b) {
        // Funzione di comparazione per ordinare i percorsi in base alla loro perdita
        return a.get_loss() < b.get_loss();
    };

    // Converti il field _path in un std::vector per usare std::sort
    vector<Path> path_vec(_path.size());
    for (unsigned int i = 0; i < _path.size(); i++) {
        path_vec[i] = _path(i);
    }

    // Ordina il vettore dei percorsi
    sort(path_vec.begin(), path_vec.end(), compare);

    // Copia il vettore ordinato di nuovo nel field _path
    for (unsigned int i = 0; i < _path.size(); i++) {
        _path(i) = path_vec[i];
    }
}

void Genetic::selection(unsigned int& parent_1, unsigned int& parent_2) {
    // Seleziona due genitori usando la selezione stocastica basata su un'esponente
    parent_1 = static_cast<unsigned int>(_population_size * pow(_rnd.Rannyu(), _sel_exp));
    parent_2 = static_cast<unsigned int>(_population_size * pow(_rnd.Rannyu(), _sel_exp));
}

void Genetic::mutation_1(vec& path) {
    // Mutazione tipo 1: scambia due città casuali nel percorso
    if (_prob_mutations(0) > _rnd.Rannyu()) {
        unsigned int index_1 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        unsigned int index_2 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        
        // Assicurati che index_1 e index_2 siano diversi
        while (index_1 == index_2) {
            index_2 = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 1));
        }

        // Scambia le città agli indici index_1 e index_2
        unsigned int temp = path(index_1);
        path(index_1) = path(index_2);
        path(index_2) = temp;

        // Verifica che il percorso sia valido
        check_path(path);
    }
}

void Genetic::mutation_2(vec& path) {
    // Mutazione tipo 2: inverti una sottosezione del percorso
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

        // Verifica che il percorso sia valido
        check_path(path);
    }
}

void Genetic::mutation_3(vec& path) {
    // Mutazione tipo 3: scambia due sotto-sezioni del percorso
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

        // Verifica che il percorso sia valido
        check_path(path);
    }
}

void Genetic::mutation_4(vec& path) {
    // Mutazione tipo 4: inverti una sezione del percorso
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

        // Verifica che il percorso sia valido
        check_path(path);
    }
}

void Genetic::crossover(vec& parent_1, vec& parent_2) {
    // Controlla se la crossover deve avvenire basandosi sulla probabilità
    if (_prob_crossover > _rnd.Rannyu()) {
        unsigned int cut_index = static_cast<unsigned int>(_rnd.Rannyu(1, _num_cities - 2));

        // Crea i figli con la parte iniziale della crossover dai genitori
        vec son_1 = parent_1.subvec(1, cut_index);
        vec son_2 = parent_2.subvec(1, cut_index);

        unsigned int counter_1 = cut_index + 1;
        unsigned int counter_2 = cut_index + 1;

        // Completa il resto di son_1 con le città mancanti da parent_2
        for (unsigned int i = 1; i < _num_cities; i++) {
            if (!any(son_1 == parent_2(i))) {
                son_1.resize(counter_1);
                son_1(counter_1 - 1) = parent_2(i);
                counter_1++;
            }
        }

        // Completa il resto di son_2 con le città mancanti da parent_1
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

        // Verifica se i percorsi dei genitori sono validi
        check_path(parent_1);
        check_path(parent_2);
    }
}

void Genetic::new_population() {
    arma::field<Path> new_gen;
    new_gen.set_size(_population_size);

    for (unsigned int i = 0; i < _half_population; i++) {
        unsigned int index_1, index_2;
        selection(index_1, index_2);
        vec parent_1 = _path(index_1).get_path();
        vec parent_2 = _path(index_2).get_path();
        crossover(parent_1, parent_2);

        // Applica mutazioni ai genitori
        mutation_1(parent_1);
        mutation_2(parent_1);
        mutation_3(parent_1);
        mutation_4(parent_1);

        mutation_1(parent_2);
        mutation_2(parent_2);
        mutation_3(parent_2);
        mutation_4(parent_2);

        // Crea nuovi percorsi per la nuova generazione
        new_gen(i).set_path_size(_num_cities + 1);
        new_gen(_population_size - 1 - i).set_path_size(_num_cities + 1);

        for (unsigned int j = 0; j < _num_cities + 1; j++) {
            new_gen(i).set_city_label(parent_1(j), j);
            new_gen(_population_size - 1 - i).set_city_label(parent_2(j), j);
        }

        // Verifica che i nuovi percorsi siano validi
        check_path(new_gen(i).get_path());
        check_path(new_gen(_population_size - 1 - i).get_path());
    }

    // Sostituisci la popolazione corrente con la nuova generazione
    for (unsigned int i = 0; i < _population_size; i++) {
        _path(i) = new_gen(i);
    }
}

void Genetic::migrate() {
    // Calcola il numero di individui da migrare
    unsigned int num_to_migrate = _population_size / 10;
    unsigned int migrants_size = num_to_migrate * (_num_cities + 1);
    vec migrants;
    migrants.set_size(migrants_size);

    // Buffer per ricevere i migranti
    vec received_migrants;
    received_migrants.set_size(migrants_size);

    // Prepara i migranti da inviare
    for (unsigned int i = 0; i < num_to_migrate; i++) {
        for (unsigned int j = 0; j < _num_cities + 1; j++) {
            migrants(i * (_num_cities + 1) + j) = this->get_path(i)(j);
        }
    }

    // Identifica i ranghi dei prossimi e precedenti processi
    int next_rank = (_rank + 1) % _size;
    int prev_rank = (_rank - 1 + _size) % _size;

    MPI_Request send_request, recv_request;
    MPI_Status send_status, recv_status;

    // Invia i migranti al prossimo processo e ricevi i migranti dal processo precedente
    MPI_Isend(migrants.memptr(), migrants_size, MPI_DOUBLE, next_rank, 0, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(received_migrants.memptr(), migrants_size, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &recv_request);

    MPI_Wait(&send_request, &send_status);
    MPI_Wait(&recv_request, &recv_status);

    // Sostituisci gli individui peggiori con quelli ricevuti
    for (unsigned int i = 0; i < num_to_migrate; i++) {
        for (unsigned int j = 0; j < _num_cities + 1; j++) {
            unsigned int label = received_migrants(i * (_num_cities + 1) + j);
            _path(_population_size - 1 - i).set_city_label(label, j);
        }
        // Verifica che il percorso ricevuto sia valido
        check_path(_path(_population_size - 1 - i).get_path());
    }

    // Ricalcola la funzione di costo per gli individui migrati
    for (unsigned int i = _population_size - 1 - num_to_migrate; i < _population_size; i++) {
        this->loss_function(i);
    }
}
