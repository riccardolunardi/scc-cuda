#include "utils.cpp"
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_FW_BW false

void f_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * colors, bool * is_visited, bool * is_eliminated, bool * is_expanded, bool &stop){
	for(int v=0; v < num_nodes; v++) {
        DEBUG_MSG("Checking " << v, "...", DEBUG_F_KERNEL);
        DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_visited[" << v << "] -> ", is_visited[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_expanded["<< v <<"] -> ", is_expanded[v], DEBUG_F_KERNEL);	
		
		if(!is_eliminated[v] && is_visited[v] && !is_expanded[v]) {
			is_expanded[v] = true;
			DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_F_KERNEL);

			for(int u = nodes[v]; u < nodes[v+1]; u++) {		
				DEBUG_MSG("		Nodo " << v << " connesso a nodo ", adjacency_list[u], DEBUG_F_KERNEL);	
				DEBUG_MSG("		is_eliminated[" << adjacency_list[u] << "] -> ", is_eliminated[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		is_visited[" << adjacency_list[u] << "] -> ", is_visited[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		colors["<<v<<"] == colors["<<adjacency_list[u]<<"] -> " << colors[v] << " == ", colors[adjacency_list[u]], DEBUG_F_KERNEL);

				if(!is_eliminated[adjacency_list[u]] && !is_visited[adjacency_list[u]] && colors[v] == colors[adjacency_list[u]]) {
					DEBUG_MSG("			is_visited[" << adjacency_list[u] << "] -> ", "TRUE", DEBUG_F_KERNEL);
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
			DEBUG_MSG("is_visited[" << i << "] -> ", is_visited[i], DEBUG_REACH);
		}
    }
}

void trimming_kernel(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * colors, bool * is_eliminated, bool &stop){
	bool elim;
	for(int v=0; v < num_nodes; v++) {
		DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
		if(!is_eliminated[v]){
			elim = true;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0 allora non va eliminato
			if(nodes[v] != nodes[v+1] && nodes_transpose[v] != nodes_transpose[v+1]){
				elim = false;
			}
			// Nel caso un arco di v faccia parte dello stesso sottografo, allora non va eliminato
			// Non serve farlo anche per la lista trasposta perchè alla fine l'if sui colors e sarebbe la stessa cosa
			DEBUG_MSG("	nodo " << v << " e nodo ", v+1, DEBUG_TRIMMING_KERNEL);
			DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_TRIMMING_KERNEL);
			for(int u = nodes[v]; u < nodes[v+1]; u++){
				DEBUG_MSG("adjacency_list[" << u << "] -> ", adjacency_list[u], DEBUG_TRIMMING_KERNEL);
				if(colors[adjacency_list[u]] == colors[v]){
					elim = false;
				}
			}
			if(elim){
				is_eliminated[v] = true;
				DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
				stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * colors, bool * is_eliminated) {
    bool stop = false;
    while(!stop) {
        stop = true;
        trimming_kernel(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, colors, is_eliminated, stop);
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
		}
    }
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
		pivots[i] = 4;
		colors[i] = 0;
	}

    bool stop = false;

    while (!stop){
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, colors, fw_is_visited, is_eliminated, fw_is_expanded);
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, colors, bw_is_visited, is_eliminated, bw_is_expanded);
		DEBUG_MSG("Trimming:" , "TESTATO SOLO TRAMITE IN/OUTDEGREE PERCHÈ SERVIREBBERO ANCHE I COLORS E NON LI ABBIAMO ANCORA FATTI", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, colors, is_eliminated);
        break;
		pivot_selection(pivots, colors, fw_is_visited, bw_is_visited, is_eliminated);
        update(colors, fw_is_visited, bw_is_visited, is_eliminated, stop);
    }

    for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", fw_is_visited[i], DEBUG_FW_BW);
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", bw_is_visited[i], DEBUG_FW_BW);
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", is_eliminated[i], DEBUG_FW_BW);
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", pivots[i], DEBUG_FW_BW);
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", colors[i], DEBUG_FW_BW);
	}
}

int main(int argc, char ** argv) {
    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	const char *filename = argv[1];
	int num_nodes, num_edges;
    int * nodes;
	int * adjacency_list;
	int * nodes_transpose;
	int * adjacency_list_transpose;

    create_graph_from_filename(filename, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);
}