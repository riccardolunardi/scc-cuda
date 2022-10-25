#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

#define DEBUG_CREATE false
#define DEBUG_MSG(str, val, print_bool){                \
    if(print_bool)                            		    \
        std::cout << str << val << std::endl;         	\
}

void read_heading_numbers(ifstream & infile, int & num_nodes, int & num_edges) {
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

void create_graph_from_header_and_stream(ifstream & infile, int num_nodes, int num_edges, int * nodes, int * adjacency_list) {
    int u, v, weight;
    string line;

    int iterator_adjacency_list = 0;

    // immagino un arco (u,v)
    // faccio la prima iterazione fuori dal while per settare il valore old_u = al primo u
    getline(infile, line);
    istringstream iss(line);
    iss >> u;
    int old_u = u;
	iss >> v;

    // debuffing, si legge finché non c'è niente
    while (iss >> weight) {}

    // dato l'arco (u,v)
    // setto il puntatore del nodo u = al posto giusto nella lista delle adiacenze
    nodes[u] = iterator_adjacency_list;
    // nella lista delle adiacenze metto v
    adjacency_list[iterator_adjacency_list] = v;

    // leggo finchè ci sono righe
    while (getline(infile, line)) {
        // Immagino un arco (u,v)
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

    // faccio puntare tutti gli ultimi nodi senza archi e il nodo dummy, ad una posizione dummy della lista d'adiacenza
    while(u < num_nodes) {
        nodes[++u] = num_edges;
    }
}

void create_transposed_graph_from_graph(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
    // scorro la lista delle adiacenze e ogni volta che trovo un nodo incremento il suo contatore
    for(int i=0; i < num_edges; i++) {
        ++ nodes_transpose[adjacency_list[i]];
    }

    int max = 0;
    // non so come cazzo spiegarlo, ma parto dal secondo e lo incremento col precedente, così facendo l'ultimo si incrementa col penultimo, che si è incrementato col second'ultimo e così via
    for(int i=0; i < num_nodes; i++) {
        max += nodes_transpose[i];
    }
    for(int i=num_nodes-1; i > -1; i--) {
        max -= nodes_transpose[i];
        nodes_transpose[i] = max;
    }

    int j=0;
    int x;

    // ora so i valori della lista nodes, ed avendo inizializzato la lista con un valore impossibile (-1) è possibile ricostruire la adjacency_list_transpose
    for(int i=0; i < num_edges; i++) {
        while(nodes[j+1] <= i && j < num_nodes - 1) {
            ++j;
        }

        x = 0;
        while(adjacency_list_transpose[nodes_transpose[adjacency_list[i]] + x] != -1) {
            ++x;
        }

        adjacency_list_transpose[nodes_transpose[adjacency_list[i]] + x] = j;
    }
}

int create_graph_from_filename(string filename, int & num_nodes, int & num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
    ifstream infile(filename);

    read_heading_numbers(infile, num_nodes, num_edges);

    DEBUG_MSG("Number of nodes: ", num_nodes, DEBUG_CREATE);
    DEBUG_MSG("Number of edges: ", num_edges, DEBUG_CREATE);

	// Definizione strutture dati principali
	nodes = new int[num_nodes];
	adjacency_list = new int[num_edges];
	nodes_transpose = new int[num_nodes];
	adjacency_list_transpose = new int[num_edges];

    // Inizializzazione delle liste 
	for (int i = 0; i < num_nodes; i++){
		nodes[i] = 0;
		nodes_transpose[i] = 0;
	}
	for (int i = 0; i < num_edges; i++){
		adjacency_list[i] = 0;
		adjacency_list_transpose[i] = -1;
	}

    create_graph_from_header_and_stream(infile, num_nodes, num_edges, nodes, adjacency_list);

    for(int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_CREATE);
    }
    for(int i = 0; i < num_edges; i++) {
        DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_CREATE);
    }

    create_transposed_graph_from_graph(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

    for(int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_CREATE);
    }
    for(int i = 0; i < num_edges; i++) {
        DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_CREATE);
    }

    return 0;
}