#include "utils.cpp"
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE true
#define DEBUG_FW_BW false
#define DEBUG_MAIN true

void f_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_visited, bool * is_eliminated, bool * is_expanded, bool &stop){
    for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_F_KERNEL);
    }

    // per ogni nodo
	for(int v=0; v < num_nodes; v++) {
        DEBUG_MSG("Checking " << v, "...", DEBUG_F_KERNEL);
        DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_visited[" << v << "] -> ", is_visited[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_expanded["<< v <<"] -> ", is_expanded[v], DEBUG_F_KERNEL);	
		
        // si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!is_eliminated[v] && is_visited[v] && !is_expanded[v]) {
            // si segna come espanso
			is_expanded[v] = true;
			DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_F_KERNEL);

            // per ogni nodo a cui punta
			for(int u = nodes[v]; u < nodes[v+1]; u++) {		
				DEBUG_MSG("		Nodo " << v << " connesso a nodo ", adjacency_list[u], DEBUG_F_KERNEL);	
				DEBUG_MSG("		is_eliminated[" << adjacency_list[u] << "] -> ", is_eliminated[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		is_visited[" << adjacency_list[u] << "] -> ", is_visited[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		pivots["<<v<<"] == pivots["<<adjacency_list[u]<<"] -> " << pivots[v] << " == ", pivots[adjacency_list[u]], DEBUG_F_KERNEL);

                // si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!is_eliminated[adjacency_list[u]] && !is_visited[adjacency_list[u]] && pivots[v] == pivots[adjacency_list[u]]) {
					DEBUG_MSG("			is_visited[" << adjacency_list[u] << "] -> ", "TRUE", DEBUG_F_KERNEL);
                    // setta il nodo puntato a visitato
					is_visited[adjacency_list[u]] = true;
                    // permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					stop = false;
				}
			}
		}
	}
}

void reach(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_visited, bool * is_eliminated, bool * is_expanded) {
    // Tutti i pivot vengono segnati come visitati
    for(int i=0; i < num_nodes; i++) {
        is_visited[ pivots[i] ] = true;
    }

    // si effettua la chiusura in avanti
    bool stop = false;
    while(!stop) {
        stop = true;
        f_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, is_visited, is_eliminated, is_expanded, stop);
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG("is_visited[" << i << "] -> ", is_visited[i], DEBUG_REACH);
		}
    }
}

void trimming_kernel(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * pivots, bool * is_eliminated, bool &stop){
	bool elim;
	for(int v=0; v < num_nodes; v++) {
		// DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
		if(!is_eliminated[v]){
			elim = true;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0 allora non va eliminato
			if(nodes[v] != nodes[v+1] && nodes_transpose[v] != nodes_transpose[v+1]){
				elim = false;
			}

            // ---- PENSO CHE QUESTO NON SIA CORRETTO (INIZIO DABRO) ----

			// Nel caso un arco di v faccia parte dello stesso sottografo, allora non va eliminato
			// Non serve farlo anche per la lista trasposta perchè alla fine l'if sui pivots e sarebbe la stessa cosa
			DEBUG_MSG("	nodo " << v << " e nodo ", v+1, DEBUG_TRIMMING_KERNEL);
			DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_TRIMMING_KERNEL);
			for(int u = nodes[v]; u < nodes[v+1]; u++){
			 	DEBUG_MSG("adjacency_list[" << u << "] -> ", adjacency_list[u], DEBUG_TRIMMING_KERNEL);
			 	// elim potrebbe ovviamente essere messo prima del for, per evitare cicli in più
				// forse va mosso, in base a come capiremo ottimizzare in CUDA
				if(pivots[adjacency_list[u]] == pivots[v] && !elim){
					DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
			 		elim = false;
				}	
			}

            // ---- PENSO CHE QUESTO NON SIA CORRETTO (FINE DABRO) ----
            // se lo è ti prego dimmi perchè dato che non ci ho capito un cazzo

			if(elim){
				is_eliminated[v] = true;
				DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
				stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * pivots, bool * is_eliminated) {
    bool stop = false;
    while(!stop) {
        stop = true;
        trimming_kernel(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, pivots, is_eliminated, stop);
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
		}
    }
}

void update(int num_nodes, int * pivots, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated, bool & stop) {
    int * write_id_for_pivots = new int[4 * num_nodes];
	for (int i = 0; i < 5 * num_nodes; i++){
		write_id_for_pivots[i] = -1;
	}

	int * colors = new int[num_nodes];
	/*
	Dai paper:
	These subgraphs are 
	1) the strongly connected component with the pivot;
	2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	3) the subgraph given by vertices in the backward closure but not in the forward closure;
	4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	they are recursively processed in parallel with the same algorithm
	*/
	stop = true;
	int new_color;
	for(int v=0; v < num_nodes; v++) {
		// Questo primo caso non ha senso di esistere, perché possiamo lasciargli il valore precedente, tanto cambiaeranno tutti gli altri
		// in realtà ha senso per conservare il valore del pivot, se poi si scopre che una volta diventato SCC il suo valore nel vettore pivots, allora il primo caso si può cancellare e moltiplicare per 3
		
		if(is_eliminated[v]){
			pivots[v] = v;
		} 
		
		if(fw_is_visited[v] == bw_is_visited[v] && fw_is_visited[v] == true){
			colors[v] = 4 * pivots[v];
		} else {
			if(fw_is_visited[v] != bw_is_visited[v] && fw_is_visited[v] == true){
				colors[v] = 4 * pivots[v] + 1;
			}else if(fw_is_visited[v] != bw_is_visited[v] && fw_is_visited[v] == false){
				colors[v] = 4 * pivots[v] + 2;
			}else if(fw_is_visited[v] == bw_is_visited[v] && fw_is_visited[v] == false){
				colors[v] = 4 * pivots[v] + 3;				
			}
				
			if(!is_eliminated[v]){
				stop = false;
				DEBUG_MSG(v, " -> non eliminato, ma non visitato da fw e bw", DEBUG_UPDATE);
			}
		}
		write_id_for_pivots[colors[v]] = v;
	}
	
	for (int i = 0; i < 4 * num_nodes; i++)
        DEBUG_MSG("write_id_for_pivots[" + to_string(i) + "] = ", write_id_for_pivots[i], DEBUG_UPDATE);
	// setto i valori dei pivot che hanno vinto la race
	// in CUDA questo è da fare dopo una sincronizzazione
	for (int i = 0; i < num_nodes; i++) {
		if(is_eliminated[i]){
			pivots[i] = i;
		}else{
			pivots[i] = write_id_for_pivots[colors[i]];
		}
	}

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_UPDATE);
}

void fw_bw(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
	bool * fw_is_visited = new bool[num_nodes];
    bool * bw_is_visited = new bool[num_nodes];
    bool * is_eliminated = new bool[num_nodes];
    bool * fw_is_expanded = new bool[num_nodes];
    bool * bw_is_expanded = new bool[num_nodes];
    int * pivots = new int[num_nodes];

	for (int i = 0; i < num_nodes; i++){
		fw_is_visited[i] = false;
		bw_is_visited[i] = false;
		is_eliminated[i] = false;
		fw_is_expanded[i] = false;
		bw_is_expanded[i] = false;
		pivots[i] = 2;
	}

    bool stop = false;

    while (!stop){
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, fw_is_visited, is_eliminated, fw_is_expanded);

        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, bw_is_visited, is_eliminated, bw_is_expanded);

		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, pivots, is_eliminated);

		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, pivots, fw_is_visited, bw_is_visited, is_eliminated, stop);
		cout << "FINE" << endl;
    }

	
    // ---- INIZIO DEBUG ----
    for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("fw_is_visited[" + to_string(i) + "] = ", fw_is_visited[i], DEBUG_FW_BW);
    for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("bw_is_visited[" + to_string(i) + "] = ", bw_is_visited[i], DEBUG_FW_BW);
    for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("is_eliminated[" + to_string(i) + "] = ", is_eliminated[i], DEBUG_FW_BW);
    for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_FW_BW);
    // ---- FINE DEBUG ----
}

int main(int argc, char ** argv) {
    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	int num_nodes, num_edges;
    int * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	bool * is_u;
    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, is_u);

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_MAIN);
	for (int i = 0; i < num_edges; i++)
        DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_MAIN);
	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_MAIN);
	for (int i = 0; i < num_edges; i++)
        DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_MAIN);
	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("is_u[" + to_string(i) + "] = ", is_u[i], DEBUG_MAIN);

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);
}