#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

#define DEBUG_UTILS true
#define DEBUG_MSG_UTILS(str, val){                     	\
    if(DEBUG_UTILS)                            		    \
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

    while (getline(infile, line)) {
        // Immagino un arco (u,v)
	    istringstream iss(line);
        iss >> u;
		iss >> v;
		//Debuffing, si legge finché non c'è niente
		while (iss >> weight) {}

        adjacency_list[++iterator_adjacency_list] = v;

        // se il nodo u è diverso dal precedente, allora aggiorna il puntatore alla lista delle adiacenze
        if(old_u != u) {
            nodes[u] = iterator_adjacency_list;
            old_u = u;
        }
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
        if (nodes_transpose[i] == 0) {
            nodes_transpose[i] = -1;
        } else {
            max -= nodes_transpose[i];
            nodes_transpose[i] = max;
        }
    }

    int i_nodes = 0;
    int i_nodes_transpose = 0;
    while(nodes[i_nodes] == -1) {++i_nodes;}
    while(nodes_transpose[i_nodes_transpose] == -1) {++i_nodes_transpose;}
    int i_adjacency_list = 0;
    int i_adjacency_list_transpose = 0;
    int find_hole_adjacency_list_traspose = 0;
    int num_minus_1 = 0;

    while(i_adjacency_list < num_edges && i_nodes < 5) {
        num_minus_1 = 0;
        while(nodes[i_nodes] < i_adjacency_list) {
            ++i_nodes;
            if(nodes[i_nodes] == -1) {
                ++num_minus_1;
            }
        }

        find_hole_adjacency_list_traspose = 0;
        while(adjacency_list_transpose[nodes_transpose[adjacency_list[i_adjacency_list]] + find_hole_adjacency_list_traspose] != -1) {
            ++find_hole_adjacency_list_traspose;
        }

        DEBUG_MSG_UTILS("-- i_nodes ", i_nodes);
        DEBUG_MSG_UTILS("------ num_minus_1 ", num_minus_1);
        adjacency_list_transpose[nodes_transpose[adjacency_list[i_adjacency_list]] + find_hole_adjacency_list_traspose] = i_nodes - num_minus_1 ;

        ++i_adjacency_list;
    }


    // ora so i valori della lista nodes_transpose, ed avendo inizializzato la lista con un valore impossibile (-1) è possibile ricostruire la adjacency_list_transpose
    // for(int i = 0; i < num_edges; i++) {
    //     while(nodes[index_nodes] <= i && index_nodes < num_nodes - 1) {
    //         ++index_nodes;
    //         DEBUG_MSG_UTILS("index_nodes = " + to_string(index_nodes) + " i = ", i);
    //     }

    //     x = 0;
    //     while(adjacency_list_transpose[nodes_transpose[adjacency_list[i]] + x] != -1) {
    //         ++x;
    //     }

    //     adjacency_list_transpose[nodes_transpose[adjacency_list[i]] + x] = index_nodes;
    // }
}

int create_graph_from_filename(string filename, int & num_nodes, int & num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
    ifstream infile(filename);

    read_heading_numbers(infile, num_nodes, num_edges);

    DEBUG_MSG_UTILS("Number of nodes: ", num_nodes);
    DEBUG_MSG_UTILS("Number of edges: ", num_edges);

	// Definizione strutture dati principali
	nodes = new int[num_nodes];
	adjacency_list = new int[num_edges];
	nodes_transpose = new int[num_nodes];
	adjacency_list_transpose = new int[num_edges];

    // Inizializzazione delle liste 
	for (int i = 0; i < num_nodes; i++){
		nodes[i] = -1;
		nodes_transpose[i] = 0;
	}
	for (int i = 0; i < num_edges; i++){
		adjacency_list[i] = 0;
		adjacency_list_transpose[i] = -1;
	}

    create_graph_from_header_and_stream(infile, num_nodes, num_edges, nodes, adjacency_list);

    for(int i = 0; i < num_nodes; i++) {
        DEBUG_MSG_UTILS("nodes[" + to_string(i) + "] = ", nodes[i]);
    }
    for(int i = 0; i < num_edges; i++) {
        DEBUG_MSG_UTILS("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i]);
    }

    create_transposed_graph_from_graph(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

    for(int i = 0; i < num_nodes; i++) {
        DEBUG_MSG_UTILS("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i]);
    }
    for(int i = 0; i < num_edges; i++) {
        DEBUG_MSG_UTILS("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i]);
    }

    return 0;
}