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

#define PRINT_RESULTS 1

void trimming_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * nodes_transpose, unsigned * adjacency_list, unsigned * adjacency_list_transpose, unsigned * pivots, char * status, bool &stop){
	// Esegue un solo ciclo di eliminazione dei nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati

	bool elim, forward;

	// Per ogni nodo v
	for(unsigned v=0; v < num_nodes; v++) {
		// Se non è stato eliminato
		if(!get_is_eliminated(status[v])){
			elim = true;
			forward = false;
			
			// Se v, contando solo i nodi non eliminati, ha sia in_degree > 0 che out_degree > 0 allora non va eliminato
			for(unsigned u = nodes[v]; u < nodes[v+1]; u++){
				if(!get_is_eliminated(status[adjacency_list[u]])) {
					forward = true;
				}
			}
			if(forward) {
				for(unsigned u = nodes_transpose[v]; u < nodes_transpose[v+1]; u++){
					if(!get_is_eliminated(status[adjacency_list_transpose[u]])) {
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

    bool stop = false;

    while(!stop) {
        stop = true;
        trimming_kernel(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status, stop);
    }
}

void reach_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status, bool &stop, bool (*get_visited)(char), bool (*get_expanded)(char), void (*set_visited)(char &), void (*set_expanded)(char &)){
	// Esecuzione di un singolo ciclo della chiusura in avanti

	// Per ogni nodo v
	for(unsigned v=0; v < num_nodes; v++) {
        // Si controlla se v non è stato eliminato, se è stato visitato e se non è stato espanso
		if(!get_is_eliminated(status[v]) && get_visited(status[v]) && !get_expanded(status[v])) {
            // Si segna come espanso, per non ricontrollarlo più nelle prossime iterazioni
			set_expanded(status[v]);

            // Per ogni nodo u a cui punta
			for(unsigned u = nodes[v]; u < nodes[v+1]; u++) {
                // Si controlla se u non è stato eliminato e se non è stato visitato
				if(!get_is_eliminated(status[adjacency_list[u]]) && !get_visited(status[adjacency_list[u]])) {
                    // Setta il nodo u come visitato
					set_visited(status[adjacency_list[u]]);
                    // Si è trovato un altro nodo visitato ancora da espandere, quindi continuo il ciclo reach
					stop = false;
				}
			}
		}
	}
}

void reach(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status, bool (*get_visited)(char), bool (*get_expanded)(char), void (*set_visited)(char &), void (*set_expanded)(char &)) {
	// Esecuzione ricorsiva della chiusura in avanti

    // Tutti i pivot vengono segnati come visitati
    for(unsigned i=0; i < num_nodes; i++) {
        set_visited(status[pivots[i]]);
    }

    // Si effettua la chiusura in avanti
    bool stop = false;
    while(!stop) {
        stop = true;
        reach_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, status, stop, get_visited, get_expanded, set_visited, set_expanded);
    }
}

void update(unsigned num_nodes, unsigned * pivots, char * status, bool & stop) {
	// Esegue l'update dei valori del pivot, facendo una race

	// Il vettore write_id_for_pivots serve per poter eleggere un pivot per ogni sotto-grafo
    unsigned * write_id_for_pivots = (unsigned*) malloc(4 * num_nodes * sizeof(unsigned));
	// Il vettore colors serve per poter distinguire i sotto-grafi
	unsigned * colors = (unsigned*) malloc(num_nodes * sizeof(unsigned));

	stop = true;

	// Per ogni nodo v
	for(unsigned v = 0; v < num_nodes; v++) {
		// Se non è stato eliminato
		if(!get_is_eliminated(status[v])){
			// Se fa parte di una SCC, quindi è stato visitato sia in avanti che all'indietro
			if(get_is_fw_visited(status[v]) == get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == true){
				colors[v] = 4 * pivots[v];
			} else {
				stop = false;

				// Se è stato visitato solo in avanti
				if(get_is_fw_visited(status[v]) != get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == true){
					colors[v] = 4 * pivots[v] + 1;
				// Se è stato visitato solo all'indietro
				}else if(get_is_fw_visited(status[v]) != get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == false){
					colors[v] = 4 * pivots[v] + 2;
				// Se non è stato visitato né in avanti né all'indietro
				}else if(get_is_fw_visited(status[v]) == get_is_bw_visited(status[v]) && get_is_fw_visited(status[v]) == false){
					colors[v] = 4 * pivots[v] + 3;	
				}
			}

			// Su questa riga viene effettuata la race
			// Ogni sotto-grafo avrà come pivot l'ultimo nodo che esegue questa riga
			write_id_for_pivots[colors[v]] = v;
		}
	}

	// Setto i valori dei pivot che hanno vinto la race
	for (unsigned v = 0; v < num_nodes; v++) {
		if(get_is_eliminated(status[v])){
			// Se sono stati eliminati e non sono una SCC, allora setta il valore del pivot uguale al nodo stesso
			if(!get_is_scc(status[v])) {
				pivots[v] = v;
			}
		}else{
			// Se non sono stati eliminati, allora setta il valore del pivot uguale al nodo che ha vinto la race
			pivots[v] = write_id_for_pivots[colors[v]];
			// I nodi che fanno parte di una SCC, vengono settati come eliminati e come SCC
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
	// Funzione che esegue l'algoritmo di Forward-Backward per la ricerca delle SCC in un grafo

    pivots = (unsigned*) malloc(num_nodes * sizeof(unsigned));

	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status);

	// Si prende come primo pivot globale il primo nodo che si riesce a trovare non eliminato 
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
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, status, get_is_fw_visited, get_is_fw_expanded, set_is_fw_visited, set_is_fw_expanded);

		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, status, get_is_bw_visited, get_is_bw_expanded, set_is_bw_visited, set_is_bw_expanded);

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, pivots, status, stop);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, status);
    }
}

void trim_u_kernel(unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status) {
	// Controlla per tutte le SCC se queste ricevono archi da nodi u
	// Se le trova setta i pivot come non facenti parte di una SCC

	for(unsigned u = 0; u < num_nodes; ++u ) {
		if(get_is_u(status[u]) == true) {
			for(unsigned v = nodes[u]; v < nodes[u+1]; ++v) {
				// Dato un arco (u,v), se u non fa parte della stessa SCC di v e u fa parte di U
				// Allora setto il pivot della SCC di v come non facente parte di una SCC
				if(pivots[u] != pivots[adjacency_list[v]]) {
					set_not_is_scc(status[pivots[adjacency_list[v]]]);
				}
			}
		}
	}
}

void trim_u_propagation(unsigned num_nodes, unsigned * pivots, char * status) {
	// Se sono presenti pivot non facenti più parte di una SCC, per la cancellazione dovuta a trim_u_kernel, 
	// propaga la cancellazione agli altri nodi della stessa SCC

	for(unsigned u = 0; u < num_nodes; ++u ) {
		if(get_is_scc(status[pivots[u]])) {
			set_is_scc(status[u]);
		} else {
			set_not_is_scc(status[u]);
		}
	}
}

void eliminate_trivial_scc(unsigned num_nodes, unsigned * pivots, char * status) {
	// Setta tutti i pivot come non facenti parte di una SCC

	for (unsigned u = 0; u < num_nodes; ++u) {
		if (pivots[u] == u) {
			set_not_is_scc(status[u]);
		}
	}
}

void trim_u(const bool profiling, unsigned num_nodes, unsigned num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * pivots, char * status, bool & is_network_valid) {
	// Elimina le SCC riceventi archi da altri nodi U non facenti parte della SCC

	// Setta i pivot delle SCC come non facenti parte di una SCC se queste ricevono archi da nodi u
	trim_u_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, status);
	// Si propaga la decisione presa dalla funzione precedente su tutti i nodi delle SCC eliminate
	trim_u_propagation(num_nodes, pivots, status);
	// Dato che l'algoritmo Forward-Backward identifica anche i singoli nodi come SCC
	// Setto tutti i pivot come non facenti parte di una SCC, facendo attenzione a non rieseguire la funzione trim_u_propagation.
	// Quindi tutte le SCC da 1 nodo saranno eliminate, mentre le SCC con 2 o più nodi verranno ancora considerate tali.
	eliminate_trivial_scc(num_nodes, pivots, status);

	if (profiling){
		is_network_valid = false;
		unsigned i = 0;
		// Al primo nodo SCC che trovo mi fermo
		while(i < num_nodes && !get_is_scc(status[i])) {
			++i;
		}
		// E se esiste una SCC si setta che la network è valida
		if(i < num_nodes) {
			is_network_valid = true;
		}
	}
}

unsigned count_distinct_scc(char status[], unsigned pivots[], unsigned n){
	// Conta quanti elementi distinti ci sono in un array

	set<unsigned> s;

	// Aggiungo un elemento al set se fa parte della SCC
	// set non permette elementi ripetuti, quindi ogni pivot comparirà una volta sola
	for(int i=0; i<n; i++) {
		if(get_is_scc(status[i])) {
        	s.insert(pivots[i]);
		}
    }
	
    return s.size();
}

void routine(const bool profiling, int num_nodes, int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
	// Funzione che printa se sono presenti SCC oppure il numero di SCC trovate

    unsigned * pivots;
	bool is_network_valid;

	// Si esegue l'algoritmo Forward-Backward per la ricerca delle SCC
	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, pivots, status);
	// Si trimmano le SCC trovate se ricevono arche dai nodi U
	trim_u(profiling, num_nodes, num_edges, nodes, adjacency_list, pivots, status, is_network_valid);

	if(profiling){
		DEBUG_MSG("", is_network_valid, PRINT_RESULTS);
	}else{
		DEBUG_MSG("Number of SCCs found: ", count_distinct_scc(status, pivots, num_nodes), PRINT_RESULTS);
	}

	free(pivots);
}
