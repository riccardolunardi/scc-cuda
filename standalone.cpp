#include "./utils/is_checked.cpp"
#include "./utils/file2graph.cpp"
#include <cstring>
#include <set>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE false
#define DEBUG_FW_BW false
#define DEBUG_MAIN false
#define DEBUG_FINAL false

#define REPEAT 1

void trimming_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * nodes_transpose, unsigned * adjacency_list, unsigned * adjacency_list_transpose, unsigned * pivots, char * status, bool &stop){
	// Esegue un solo ciclo di eliminazione dei nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			status	=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati
	// @return:	status	=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati, aggiornata dopo l'esecuzione di trimming_kernel

	bool elim, forward, backward;
	for(unsigned v=0; v < num_nodes; v++) {
		if(!get_is_eliminated(status[v])){
			elim = true;
			forward = false;
			backward = false;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0, tra i soli nodi non eliminati, allora non va eliminato
			for(unsigned u = nodes[v]; u < nodes[v+1]; u++){
				if(!get_is_eliminated(status[adjacency_list[u]])) {
				// if(!get_is_eliminated(status[adjacency_list[u]]) && pivots[v] == pivots[u]) {
					forward = true;
				}
			}
			if(forward) {
				for(unsigned u = nodes_transpose[v]; u < nodes_transpose[v+1]; u++){
					if(!get_is_eliminated(status[adjacency_list_transpose[u]])) {
					// if(!get_is_eliminated(status[adjacency_list_transpose[u]]) && pivots[v] == pivots[u]) {
						elim = false;
					}
				}
			}

			if(elim){
				set_is_eliminated(status[v]);
				stop = false;
			}
		}
	}
}

void trimming(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * nodes_transpose, unsigned * adjacency_list, unsigned * adjacency_list_transpose, unsigned * pivots, char * status) {
	// Elimina ricorsivamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati
	// @return:	status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati

    bool stop = false;

    while(!stop) {
        stop = true;
        trimming_kernel(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status, stop);
    }
}

void reach_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status, bool &stop, bool (*get_visited)(char), bool (*get_expanded)(char), void (*set_visited)(char &), void (*set_expanded)(char &)){
	// Esecuzione di un singolo ciclo della chiusura in avanti
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati
	// @return 	status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati, aggiornata dopo l'esecuzione del trimming

    // for (unsigned i = 0; i < num_nodes; i++){
    //     DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_F_KERNEL);
    // }

	// Per ogni nodo
	for(unsigned v=0; v < num_nodes; v++) {
        // DEBUG_MSG("Checking " << v, "...", DEBUG_F_KERNEL);
		
        // Si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!get_is_eliminated(status[v]) && get_visited(status[v]) && !get_expanded(status[v])) {
            // Si segna come espanso
			set_expanded(status[v]);
			// DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_F_KERNEL);

            // Per ogni nodo a cui punta
			for(unsigned u = nodes[v]; u < nodes[v+1]; u++) {		
				// DEBUG_MSG("		Nodo " << v << " connesso a nodo ", adjacency_list[u], DEBUG_F_KERNEL);
				// DEBUG_MSG("		pivots["<<v<<"] == pivots["<<adjacency_list[u]<<"] -> " << pivots[v] << " == ", pivots[adjacency_list[u]], DEBUG_F_KERNEL);

                // Si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!get_is_eliminated(status[adjacency_list[u]]) && !get_visited(status[adjacency_list[u]])) {
				// if(!get_is_eliminated(status[adjacency_list[u]]) && !get_visited(status[adjacency_list[u]]) && pivots[v] == pivots[adjacency_list[u]]) {
					// DEBUG_MSG("			status[" << adjacency_list[u] << "] -> ", "TRUE", DEBUG_F_KERNEL);
                    // Setta il nodo puntato a visitato
					set_visited(status[adjacency_list[u]]);
                    // Permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					stop = false;
				}
			}
		}
	}
}

void reach(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status, bool (*get_visited)(char), bool (*get_expanded)(char), void (*set_visited)(char &), void (*set_expanded)(char &)) {
	// Esecuzione ricorsiva della chiusura in avanti
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati
	// @return 	status		=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati, aggiornata dopo l'esecuzione del trimming

    // Tutti i pivot vengono segnati come visitati
    for(unsigned i=0; i < num_nodes; i++) {
        set_visited(status[pivots[i]]);
    }

    // si effettua la chiusura in avanti
    bool stop = false;
    while(!stop) {
        stop = true;
        reach_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, status, stop, get_visited, get_expanded, set_visited, set_expanded);
    }
}

void update(unsigned num_nodes, unsigned * pivots, char * status, bool & stop) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			status		= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return: pivots			= Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC, aggiornata dopo l'esecuzione di update

    unsigned * write_id_for_pivots = (unsigned*) malloc(4 * num_nodes * sizeof(unsigned));
	unsigned * colors = (unsigned*) malloc(num_nodes * sizeof(unsigned));

	// Dai paper:
	// These subgraphs are 
	// 		1) the strongly connected component with the pivot;
	// 		2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	// 		3) the subgraph given by vertices in the backward closure but not in the forward closure;
	// 		4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	// The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	// they are recursively processed in parallel with the same algorithm

	stop = true;
	for(unsigned v = 0; v < num_nodes; v++) {
		if(!get_is_eliminated(status[v])){
			if(get_is_fw_visited(status[v]) == get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == true){
				colors[v] = 4 * pivots[v];
			} else {
				stop = false;

				if(get_is_fw_visited(status[v]) != get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == true){
					colors[v] = 4 * pivots[v] + 1;
				}else if(get_is_fw_visited(status[v]) != get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == false){
					colors[v] = 4 * pivots[v] + 2;
				}else if(get_is_fw_visited(status[v]) == get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == false){
					colors[v] = 4 * pivots[v] + 3;	
				}
			}

			write_id_for_pivots[colors[v]] = v;
		}
	}

	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	for (unsigned v = 0; v < num_nodes; v++) {
		if(get_is_eliminated(status[v])){
			if(!get_is_scc(status[v])) {
				pivots[v] = v;
			}
		}else{
			pivots[v] = write_id_for_pivots[colors[v]];
			if(colors[v] % 4 == 0) {
				set_is_eliminated(status[v]);
				set_is_scc(status[v]);
			}
		}
	}

	free(colors);
	free(write_id_for_pivots);
}

void fw_bw(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, unsigned *& pivots, char * status) {
	// Calcola il Forward-Backward di un grafo

    pivots = (unsigned*) malloc(num_nodes * sizeof(unsigned));
	char * bw_status = (char*) malloc(num_nodes * sizeof(char));

	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status);

	// Si prende come primo pivot globale, il primo nodo che si riesce a trovare non eliminato 
	unsigned v = 0;
	while(v < num_nodes && get_is_eliminated(status[v])) {
		++v;
	}
	for (unsigned i = 0; i < num_nodes; i++){
		pivots[i] = v;
	}

	DEBUG_MSG("pivot = " , v, DEBUG_FW_BW);

    bool stop = false;

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
    while (!stop){
		memcpy(bw_status, status, num_nodes);

		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, status, get_is_fw_visited, get_is_fw_expanded, set_is_fw_visited, set_is_fw_expanded);

		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, bw_status, get_is_bw_visited, get_is_bw_expanded, set_is_bw_visited, set_is_bw_expanded);

		for(int i=0; i<num_nodes; i++) {
			status[i] |= bw_status[i];
		}

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, pivots, status, stop);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status);
    }
}

void trim_u_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status) {
	// Setta i pivot delle SCC uguale a -1 se questi ricevono archi da nodi u
	// param: 	pivots = 	Lista che per ogni 'v' dice il valore del pivot della SCC

	for(unsigned u = 0; u < num_nodes; ++u ) {
		if(get_is_u(status[u]) == true) {
			for(unsigned v = nodes[u]; v < nodes[u+1]; ++v) {
				if(pivots[u] != pivots[adjacency_list[v]]) {
					set_not_is_scc(status[pivots[adjacency_list[v]]]);
				}
			}
		}
	}
}

void trim_u_propagation(unsigned num_nodes, unsigned * pivots, char * status) {
	// Se alcuni pivot sono settati a -1, per la cancellazione dovuta a collegamenti con nodi u, 
	// propaga la cancellazione agli altri membri della SCC
	// param: 	pivots = 	Lista contenente i pivot delle SCC

	for(unsigned u = 0; u < num_nodes; ++u ) {
		if(get_is_scc(status[pivots[u]])) {
			set_is_scc(status[u]);
		} else {
			set_not_is_scc(status[u]);
		}
	}
}

void eliminate_trivial_scc(unsigned num_nodes, unsigned * pivots, char * status) {
	for (unsigned u = 0; u < num_nodes; ++u) {
		if (pivots[u] == u) {
			set_not_is_scc(status[u]);
		}
	}
}

void is_scc_adjust_host(unsigned num_nodes, unsigned * pivots, char * status) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// In questa fase la lista ha -1 nei valori dei pivot. Per fixare, i nodi facendi parte di quella SCC
	// andranno a scrivere nella posizione del pivot, il valore del pivot stesso

	for (unsigned u = 0; u < num_nodes; ++u) {
		if (!get_is_scc(status[u])) {
			set_not_is_scc(status[pivots[u]]);
		}
	}
}

void trim_u(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status) {
	// Elimina le SCC riceventi archi da altri nodi U non facenti parte della SCC
	// @param:	pivots 	=	Lista che per ogni 'v' dice il valore del pivot della SCC
	// 			status	=	Lista che per ogni 'v' contiene 8 bit che rappresentano degli stati

	trim_u_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, status);
	trim_u_propagation(num_nodes, pivots, status);
	eliminate_trivial_scc(num_nodes, pivots, status);

	// Si prende come pivot, il primo pivot che si riesce a trovare facente parte di una scc 
	unsigned pivot_riferimento;
	bool pivot_riferimento_found = false;
	unsigned v = 0;
	while(v < num_nodes && !get_is_scc(status[v])) {
		++v;
	}

	// Se trova una SCC, prende il pivot del primo nodo che fa parte di essa come pivot riferimento
	// altrimenti setta il pivot di riferimento a -1 (impossibile)
	if(v < num_nodes) {
		pivot_riferimento = pivots[v];
		pivot_riferimento_found = true;
	}

	// salviamo in una struttura dati tutti i nodi aventi il pivot_riferimento
	if(pivot_riferimento_found){
		set<unsigned> nodi_in_scc;
		for (unsigned u = 0; u < num_nodes; ++u) {
			if (pivots[u] == pivot_riferimento) {
				nodi_in_scc.insert(u);
			}
		}

		// stampa ogni elemento di nodi_in_scc
		for (auto it = nodi_in_scc.begin(); it != nodi_in_scc.end(); ++it) {
			cout << *it << " ";
		}
	}else{
		cout << "No SCCs found";
	}

	cout << endl;
}

unsigned count_distinct_scc(char status[], unsigned pivots[], unsigned n){
	// Conta quanti elementi distinti ci sono in un array
	// @param:	arr =	Array in cui contare il numero di elementi diverso
	// 			n 	=	Numero di elementi nell'array
	// @return:	res =	Numero di elementi diversi nell'array

	set<unsigned> s;

	for(int i=0; i<n; i++) {
		if(get_is_scc(status[i])) {
        	s.insert(pivots[i]);
		}
    }
	
    return s.size();
}

int routine(int num_nodes, int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
    unsigned * pivots;

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, pivots, status);
	trim_u(num_nodes, num_edges, nodes, adjacency_list, pivots, status);

	free(pivots);
	return 0;
}

int main(int argc, char ** argv) {
    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	unsigned num_nodes, num_edges;
    unsigned * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	char * status;

    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);

	char * og_status;
	og_status = (char *) malloc(num_nodes * sizeof(char));
	memcpy(og_status, status, num_nodes);

	for(int i=0;i<REPEAT;i++){
		memcpy(status, og_status, num_nodes);
		routine(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

	free(nodes);
	free(adjacency_list);
	free(nodes_transpose);
	free(adjacency_list_transpose);
	free(status);
	free(og_status);
}