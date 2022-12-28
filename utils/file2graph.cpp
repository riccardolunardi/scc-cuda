#ifndef FILE2GRAPH
#define FILE2GRAPH

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <cstring>
#include <numeric>
#include "../utils/is_checked.cpp"
using namespace std;

#define DEBUG_CREATE false
#define DEBUG_MSG(str, val, print_bool){                \
    if(print_bool)                            		    \
        std::cout << str << val << std::endl;         	\
}

void read_heading_numbers(ifstream & infile, unsigned & num_nodes, unsigned & num_edges) {
	string line;
    char percentage_sign;
	getline(infile, line);
	istringstream iss(line);

    // scarto il carattere % che serve a capire quale è la riga d'intestazione
	iss >> percentage_sign;
	iss >> num_edges;
	iss >> num_nodes;

    // aggiungo il nodo dummy
    ++num_nodes;
}

void create_graph_from_header_and_stream(ifstream & infile, unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, char * status) {
    unsigned u, v, weight;
    string line;

    unsigned iterator_adjacency_list = 0;

    // immagino un arco (u,v)
    // faccio la prima iterazione fuori dal while per settare il valore old_u = al primo u
    getline(infile, line);
    istringstream iss(line);
    iss >> u;
    unsigned old_u = u;
	iss >> v;

    // debuffing, si legge finché non c'è niente
    while (iss >> weight) {}

    // dato l'arco (u,v)
    // setto il puntatore del nodo u = al posto giusto nella lista delle adiacenze
    nodes[u] = iterator_adjacency_list;
    // nella lista delle adiacenze metto v
    adjacency_list[iterator_adjacency_list] = v;

    // Vecchio codice meno performante
    /*for(unsigned i = 0; i < num_edges - 1; ++i) {
        // Immagino un arco (u,v)
        getline(infile, line);
	    istringstream iss(line);
        iss >> u;
		iss >> v;

		// Debuffing, si legge finché non c'è niente
		while (iss >> weight) {}

        adjacency_list[++iterator_adjacency_list] = v;

        // se il nodo u è diverso dal precedente, setto tutti i valori di nodes da [old_u + 1] a [u] compreso con il valore dell'iteratore della lista d'adiacenza
        if(old_u != u) {
            while(old_u < u) {
                nodes[++old_u] = iterator_adjacency_list;
            }
        }
    }
    */

    for(unsigned i = 0; i < num_edges - 1; ++i) {
        infile >> u >> v;
        adjacency_list[++iterator_adjacency_list] = v;

        if (old_u != u) {
            while (old_u < u) {
                nodes[++old_u] = iterator_adjacency_list;
            }
        }
    }

    // faccio puntare tutti gli ultimi nodi senza archi e il nodo dummy, ad una posizione dummy della lista d'adiacenza
    while(u < num_nodes - 1) {
        nodes[++u] = num_edges;
    }

    // leggo finchè ci sono righe
    while (infile >> u) {
        status[u] = (status[u] | 32) & 251;
        
        //set_not_is_eliminated(status[u]);
        //set_is_u(status[u]);
    }
}

void create_transposed_graph_from_graph(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose) {
    // creo una lista per capire se un nodo nella trapsose adjacency list è occupato
    // All'inizio, setto tutti a non occupati per iniziare
	bool * adjacency_list_transpose_occupied = (bool*) malloc(num_edges * sizeof(bool));
    memset(adjacency_list_transpose_occupied, 0, num_edges);

    // scorro la lista delle adiacenze e ogni volta che trovo un nodo incremento il suo contatore
    for(unsigned i=0; i < num_edges; i++) {
        ++ nodes_transpose[adjacency_list[i]];
    }

    // Prima lo faceva un for, ma così è più performante
    unsigned max = std::accumulate(nodes_transpose, nodes_transpose + num_nodes, 0);

    // Parto dall'ultimo nodo e computo la posizione inizale nella lista di adiacenza, procedo a ritroso nel vettore
    for(unsigned i = num_nodes; i > 0; i--) {
        max -= nodes_transpose[i-1];
        nodes_transpose[i-1] = max;
    }

    unsigned pointed_node = 0;
    unsigned first_position_available;

    // ora so i valori della lista nodes_transpose, ed avendo inizializzato la adjacency_list_transpose_occupied con valore false so dove posizionare i valori nella adjacency_list_transpose
    for(unsigned index_adjacency_list = 0; index_adjacency_list < num_edges; ++index_adjacency_list) {
        // cerco i nodi che vengono puntati da altri nodi, quando esce dal ciclo pointed_node avrà il valore di un nodo puntato dal nodo nodes[pointed_node]
        while(nodes[pointed_node+1] <= index_adjacency_list && pointed_node < num_nodes - 1) {
            ++pointed_node;
        }

        // trovo il primo posto disponibile nella corretta posizione della adjacency_list_transpose
        first_position_available = 0;
        // se il nodo nella lista trasposta delle adiacenze è occupato, passo al nodo successivo
        while(adjacency_list_transpose_occupied[nodes_transpose[adjacency_list[index_adjacency_list]] + first_position_available]) {
            ++first_position_available;
        }

        // metto nella prima posizione libera, del corretto nodo, della adjacency_list_transpose il nodo puntato nel grafo trasposto
        adjacency_list_transpose[nodes_transpose[adjacency_list[index_adjacency_list]] + first_position_available] = pointed_node;
        adjacency_list_transpose_occupied[nodes_transpose[adjacency_list[index_adjacency_list]] + first_position_available] = true;
    }

    free(adjacency_list_transpose_occupied);
}

void create_graph_from_filename(string filename, unsigned & num_nodes, unsigned & num_edges, unsigned *& nodes, unsigned *& adjacency_list, unsigned *& nodes_transpose, unsigned *& adjacency_list_transpose, char *& status) {
    ifstream infile(filename);

    read_heading_numbers(infile, num_nodes, num_edges);

	// Definizione strutture dati principali
	nodes = (unsigned*) malloc(num_nodes * sizeof(unsigned));
	adjacency_list = (unsigned*) malloc(num_edges * sizeof(unsigned));
	nodes_transpose = (unsigned*) malloc(num_nodes * sizeof(unsigned));
	adjacency_list_transpose = (unsigned*) malloc(num_edges * sizeof(unsigned));
    status = (char *) malloc(num_nodes * sizeof(char));
    // Li setto tutti a eliminati
    // Quando troverò i nodi di U, setterò gli stessi come non eliminati
	memset(status, 4, num_nodes); // Il vecchio valore era 68

    // Inizializzazione delle liste 
    memset(nodes, 0, num_nodes * sizeof(unsigned));
    memset(adjacency_list, 0, num_edges * sizeof(unsigned));
    memset(nodes_transpose, 0, num_nodes * sizeof(unsigned));

    create_graph_from_header_and_stream(infile, num_nodes, num_edges, nodes, adjacency_list, status);

    create_transposed_graph_from_graph(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

    // si vuole eliminare dalla vista il nodo dummy, tenendo in considerazione che logicamente è presente
    --num_nodes;

    DEBUG_MSG("Number of nodes: ", num_nodes, DEBUG_CREATE);
    DEBUG_MSG("Number of edges: ", num_edges, DEBUG_CREATE);

    if (DEBUG_CREATE){
        for(int i = 0; i < num_nodes; i++)
            DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_CREATE);
        for(int i = 0; i < num_edges; i++)
            DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_CREATE);
        for(int i = 0; i < num_nodes; i++)
            DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_CREATE);
        for(int i = 0; i < num_edges; i++)
            DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_CREATE);
    }
}

#endif