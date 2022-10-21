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

void create_graph_from_file(ifstream & infile, int num_nodes, int num_edges, int * nodes, int * adjacency_list) {
    int u, v, weight;
    string line;

    int i = 0;
    int j = -1;

    int old_u = 0;
    while (std::getline(infile, line)) {
        // Immagino un arco (u,v)
	    istringstream iss(line);
        iss >> u;
		iss >> v;
		//Debuffing, si legge finché non c'è niente
		while (iss >> weight) {}

        adjacency_list[++j] = v;

        // se il nodo u è diverso dal precedente, allora aggiorna il puntatore alla lista delle adiacenze
        if(old_u != u) {
            nodes[u] = j;
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

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}
	const char *filename = argv[1];
    ifstream infile(filename);
	int num_nodes, num_edges;

    read_heading_numbers(infile, num_nodes, num_edges);

    DEBUG_MSG_UTILS("Number of nodes: ", num_nodes);
    DEBUG_MSG_UTILS("Number of edges: ", num_edges);

	// Definizione strutture dati principali
	int *nodes = new int[num_nodes];
	int *adjacency_list = new int[num_edges];
	int *nodes_transpose = new int[num_nodes];
	int *adjacency_list_transpose = new int[num_edges];

    // Inizializzazione delle liste 
	for (int i = 0; i < num_nodes; i++){
		nodes[i] = 0;
		nodes_transpose[i] = 0;
	}
	for (int i = 0; i < num_edges; i++){
		adjacency_list[i] = 0;
		adjacency_list_transpose[i] = -1;
	}

    create_graph_from_file(infile, num_nodes, num_edges, nodes, adjacency_list);

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
}