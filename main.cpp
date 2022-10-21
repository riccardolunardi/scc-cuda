#include "adj_list.cpp"
#include "reverse_adj_list.cpp"

#define debug 1
using namespace std;

#ifdef debug
#define DEBUG_MSG(str, value)                 \
	do                                 \
	{                                  \
		std::cout << str << value << std::endl; \
	} while (false)
#else
#define DEBUG_MSG(str, value) \
	do                 \
	{                  \
	} while (false)
#endif

void reach(int * num_nodes, int * num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_visited) {

}

void trimming(int * num_nodes, int * num_edges, int * nodes, int * adjacency_list, bool * is_eliminated) {

}

void pivot_selection(int * pivots, int * colors, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated) {

}

void update(int * colors, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated, bool & stop) {

}

void fw_bw(int * num_nodes, int * num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
	bool * fw_is_visited = new bool[num_nodes];
    bool * bw_is_visited = new bool[num_nodes];
    bool * is_eliminated = new bool[num_nodes];
    int * pivots = new int[num_nodes];
    int * colors = new int[num_nodes];

	for (int i = 0; i < num_nodes; i++){
		fw_is_visited[i] = false;
		bw_is_visited[i] = false;
		is_eliminated[i] = false;
		pivots[i] = 0;
		colors[i] = i;
	}

    bool stop = false;

    // while (!stop){
    //     reach(num_nodes, num_edges, nodes, adjacency_list, pivots, fw_is_visited);
    //     reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, bw_is_visited);
    //     trimming(num_nodes, num_edges, nodes, adjacency_list, is_eliminated);
    //     pivot_selection(pivots, colors, fw_is_visited, bw_is_visited, is_eliminated);
    //     update(colors, fw_is_visited, bw_is_visited, is_eliminated, stop);
    // }

    for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG("fw_is_visited[%d]: %d", i, fw_is_visited[i]);
        DEBUG_MSG("bw_is_visited[%d]: %d", i, bw_is_visited[i]);
        DEBUG_MSG("is_eliminated[%d]: %d", i, is_eliminated[i]);
        DEBUG_MSG("pivots[%d]: %d", i, pivots[i]);
        DEBUG_MSG("colors[%d]: %d", i, colors[i]);
	}
}


int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}
	const char *filename = argv[1];
	char percentage_sign;
	int num_edges, num_nodes;

    /* ---- INIZIO LETTURA ---- */
	//Così sembra che la prima riga sia letteralmente scartata
	std::string line;
	std::ifstream infile(filename);
	std::getline(infile, line);
	std::getline(infile, line);
	std::istringstream iss(line);
    //Questo dovrebbe essere il simbolo "%" nei file, che viene acquisito qui, ma poi mai più usato
	iss >> percentage_sign;
	iss >> num_edges;
	iss >> num_nodes;
	infile.close();

    DEBUG_MSG("Number of nodes: ", num_nodes);
    DEBUG_MSG("Number of edges: ", num_edges);
    /* ---- FINE LETTURA ---- */

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
		adjacency_list_transpose[i] = 0;
	}

	// Creazione delle liste di adiacenza
	adj_list(filename, nodes, adjacency_list);

	// Creazione delle liste di adiacenza del grafo trasposto (per la backward clousure)
	// Forse si può evitare la ripetizione di codice usando il codice di adj_list leggermente modificato
	reverse_adj_list(filename, nodes_transpose, adjacency_list_transpose);

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);
}