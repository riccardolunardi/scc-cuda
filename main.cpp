#include "utils.cpp"

#define debug_main true
using namespace std;

#define DEBUG_MSG_MAIN(str, value){                     \
    if(debug_main)                                      \
        std::cout << str << value << std::endl;         \
}

void f_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * colors, bool * is_visited, bool * is_eliminated, bool * is_expanded, bool &stop){
	int read_until = 0;
	for(int v=0; v < num_nodes; v++) {
		if(debug_main){
			DEBUG_MSG_MAIN("Checking " << v, "...");
			DEBUG_MSG_MAIN("is_eliminated[" << v << "] -> ", is_eliminated[v]);
			DEBUG_MSG_MAIN("is_visited[" << v << "] -> ", is_visited[v]);
			DEBUG_MSG_MAIN("is_expanded["<< v <<"] -> ", is_expanded[v]);	
		}
		
		if(!is_eliminated[v] && is_visited[v] && !is_expanded[v]) {
			is_expanded[v] = true;

			if(v < num_nodes - 1){
				read_until = nodes[v+1];
			}else{
				read_until = num_edges;
			}
			DEBUG_MSG_MAIN("	u va da " << nodes[v] << " a ", read_until - 1);

			for(int u = nodes[v]; u < read_until; u++) {		
				if (debug_main){
				DEBUG_MSG_MAIN("		Nodo " << v << " connesso a nodo ", adjacency_list[u]);	
				DEBUG_MSG_MAIN("		is_eliminated[" << adjacency_list[u] << "] -> ", is_eliminated[adjacency_list[u]]);
				DEBUG_MSG_MAIN("		is_visited[" << adjacency_list[u] << "] -> ", is_visited[adjacency_list[u]]);
				DEBUG_MSG_MAIN("		colors["<<v<<"] == colors["<<adjacency_list[u]<<"] -> " << colors[v] << " == ", colors[adjacency_list[u]]);
				}
				if(!is_eliminated[adjacency_list[u]] && !is_visited[adjacency_list[u]] && colors[v] == colors[adjacency_list[u]]) {
					DEBUG_MSG_MAIN("			is_visited[" << adjacency_list[u] << "] -> ", "TRUE");
					is_visited[adjacency_list[u]] = true;
					stop = false;
				}
			}
		}
	}
}

void reach(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, int * colors, bool * is_visited, bool * is_eliminated, bool * is_expanded) {
    // Tutti i pivot vengono segnati come visitati
    for(int i=0; i < num_nodes; i++) {
        is_visited[ pivots[i] ] = true;
    }

    // si effettua la chiusura in avanti
    bool stop = false;
    while(!stop) {
        stop = true;
        f_kernel(num_nodes, num_edges, nodes, adjacency_list, colors, is_visited, is_eliminated, is_expanded, stop);
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG_MAIN("is_visited[" << i << "] -> ", is_visited[i]);
		}
    }
}

void trimming(int num_nodes, int num_edges, int * nodes, int * adjacency_list, bool * is_eliminated) {

}

void pivot_selection(int * pivots, int * colors, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated) {

}

void update(int * colors, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated, bool & stop) {

}

void fw_bw(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
	bool * fw_is_visited = new bool[num_nodes];
    bool * bw_is_visited = new bool[num_nodes];
    bool * is_eliminated = new bool[num_nodes];
    bool * fw_is_expanded = new bool[num_nodes];
    bool * bw_is_expanded = new bool[num_nodes];
    int * pivots = new int[num_nodes];
    int * colors = new int[num_nodes];

	for (int i = 0; i < num_nodes; i++){
		fw_is_visited[i] = false;
		bw_is_visited[i] = false;
		is_eliminated[i] = false;
		fw_is_expanded[i] = false;
		bw_is_expanded[i] = false;
		pivots[i] = 7;
		colors[i] = 0;
	}

    bool stop = false;

    while (!stop){
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, colors, fw_is_visited, is_eliminated, fw_is_expanded);
		break;
        reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, colors, bw_is_visited, is_eliminated, bw_is_expanded);
        trimming(num_nodes, num_edges, nodes, adjacency_list, is_eliminated);
        pivot_selection(pivots, colors, fw_is_visited, bw_is_visited, is_eliminated);
        update(colors, fw_is_visited, bw_is_visited, is_eliminated, stop);
    }

    for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG_MAIN("fw_is_visited[" + to_string(i) +"]", to_string(fw_is_visited[i]));
        DEBUG_MSG_MAIN("bw_is_visited[" + to_string(i) +"]", to_string(bw_is_visited[i]));
        DEBUG_MSG_MAIN("is_eliminated[" + to_string(i) +"]", to_string(is_eliminated[i]));
        DEBUG_MSG_MAIN("pivots[" + to_string(i) +"]", to_string(pivots[i]));
        DEBUG_MSG_MAIN("colors[" + to_string(i) +"]", to_string(colors[i]));
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

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);
}