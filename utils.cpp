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

void create_graph_from_header_and_stream(ifstream & infile, int num_nodes, int num_edges, int *& nodes, int *& adjacency_list, bool *& is_eliminated) {
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

    for(int i = 0; i < num_edges - 1; ++i) {
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

    // faccio puntare tutti gli ultimi nodi senza archi e il nodo dummy, ad una posizione dummy della lista d'adiacenza
    while(u < num_nodes) {
        nodes[++u] = num_edges;
    }

    // leggo finchè ci sono righe
    while (getline(infile, line)) {
	    istringstream iss(line);
        iss >> u;
        is_eliminated[u] = false;
    }
}

void create_transposed_graph_from_graph(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int *& nodes_transpose, int *& adjacency_list_transpose) {
    // scorro la lista delle adiacenze e ogni volta che trovo un nodo incremento il suo contatore
    for(int i=0; i < num_edges; i++) {
        ++ nodes_transpose[adjacency_list[i]];
    }

    int max = 0;
    // faccio la somma di tutti i valori contenuti in nodes
    for(int i=0; i < num_nodes; i++) {
        max += nodes_transpose[i];
    }
    // Parto dall'ultimo nodo e computo la posizione inizale nella lista di adiacenza, procedo a ritroso nel vettore
    for(int i=num_nodes-1; i > -1; i--) {
        max -= nodes_transpose[i];
        nodes_transpose[i] = max;
    }

    int pointed_node = 0;
    int first_position_available;

    // ora so i valori della lista nodes_transpose, ed avendo inizializzato la adjacency_list_transpose con un valore impossibile (-1) è possibile ricostruirla
    for(int index_adjacency_list = 0; index_adjacency_list < num_edges; ++index_adjacency_list) {
        // cerco i nodi che vengono puntati da altri nodi, quando esce dal ciclo pointed_node avrà il valore di un nodo puntato dal nodo nodes[pointed_node]
        while(nodes[pointed_node+1] <= index_adjacency_list && pointed_node < num_nodes - 1) {
            ++pointed_node;
        }

        // trovo il primo posto disponibile nella corretta posizione della adjacency_list_transpose
        first_position_available = 0;
        while(adjacency_list_transpose[nodes_transpose[adjacency_list[index_adjacency_list]] + first_position_available] != -1) {
            ++first_position_available;
        }

        // metto nella prima posizione libera, del corretto nodo, della adjacency_list_transpose il nodo puntato nel grafo trasposto
        adjacency_list_transpose[nodes_transpose[adjacency_list[index_adjacency_list]] + first_position_available] = pointed_node;
    }
}

int create_graph_from_filename(string filename, int & num_nodes, int & num_edges, int *& nodes, int *& adjacency_list, int *& nodes_transpose, int *& adjacency_list_transpose, bool *& is_eliminated) {
    ifstream infile(filename);

    read_heading_numbers(infile, num_nodes, num_edges);

	// Definizione strutture dati principali
	nodes = new int[num_nodes];
	adjacency_list = new int[num_edges];
	nodes_transpose = new int[num_nodes];
	adjacency_list_transpose = new int[num_edges];
	is_eliminated = new bool[num_nodes];

    // Inizializzazione delle liste 
	for (int i = 0; i < num_nodes; i++){
		nodes[i] = 0;
		nodes_transpose[i] = 0;
	    is_eliminated[i] = true;
	}
	for (int i = 0; i < num_edges; i++){
		adjacency_list[i] = 0;
		adjacency_list_transpose[i] = -1;
	}

    create_graph_from_header_and_stream(infile, num_nodes, num_edges, nodes, adjacency_list, is_eliminated);

    create_transposed_graph_from_graph(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

    // si vuole eliminare dalla vista il nodo dummy, tenendo in considerazione che logicamente è presente
    --num_nodes;

    DEBUG_MSG("Number of nodes: ", num_nodes, DEBUG_CREATE);
    DEBUG_MSG("Number of edges: ", num_edges, DEBUG_CREATE);

    for(int i = 0; i < num_nodes; i++)
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_CREATE);
    for(int i = 0; i < num_edges; i++)
        DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_CREATE);
    for(int i = 0; i < num_nodes; i++)
        DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_CREATE);
    for(int i = 0; i < num_edges; i++)
        DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_CREATE);

    return 0;
}