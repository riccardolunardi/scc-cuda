#include "utils.cpp"
#include <cstring>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE false
#define DEBUG_FW_BW false
#define DEBUG_MAIN false
#define DEBUG_FINAL false

void trimming_kernel(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * adjacency_list_transpose, int * pivots, bool * is_eliminated, bool &stop){
	// Esegue un solo ciclo di eliminazione dei nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming

	bool elim, forward, backward;
	for(int v=0; v < num_nodes; v++) {
		DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
		if(!is_eliminated[v]){
			elim = true;
			forward = false;
			backward = false;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0, tra i soli nodi non eliminati, allora non va eliminato
			for(int u = nodes[v]; u < nodes[v+1]; u++){
				if(!is_eliminated[adjacency_list[u]]) {
					forward = true;
				}
			}
			if(forward) {
				for(int u = nodes_transpose[v]; u < nodes_transpose[v+1]; u++){
					if(!is_eliminated[adjacency_list_transpose[u]]) {
						backward = true;
					}
				}
			}
			if(backward) {
				elim = false;
			}

			if(elim){
				is_eliminated[v] = true;
				DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_TRIMMING_KERNEL);
				stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int num_edges, int * nodes, int * nodes_transpose, int * adjacency_list, int * adjacency_list_transpose, int * pivots, bool * is_eliminated) {
	// Elimina ricorsivamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming

    bool stop = false;
	for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	}
    while(!stop) {
        stop = true;
        trimming_kernel(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, is_eliminated, stop);
    }

	for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	}
}

void reach_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_visited, bool * is_eliminated, bool * is_expanded, bool &stop){
	// Esecuzione di un singolo ciclo della chiusura in avanti
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming

    for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_F_KERNEL);
    }

	// Per ogni nodo
	for(int v=0; v < num_nodes; v++) {
        DEBUG_MSG("Checking " << v, "...", DEBUG_F_KERNEL);
        DEBUG_MSG("is_eliminated[" << v << "] -> ", is_eliminated[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_visited[" << v << "] -> ", is_visited[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_expanded["<< v <<"] -> ", is_expanded[v], DEBUG_F_KERNEL);	
		
        // Si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!is_eliminated[v] && is_visited[v] && !is_expanded[v]) {
            // Si segna come espanso
			is_expanded[v] = true;
			DEBUG_MSG("	u va da " << nodes[v] << " a ", nodes[v+1], DEBUG_F_KERNEL);

            // Per ogni nodo a cui punta
			for(int u = nodes[v]; u < nodes[v+1]; u++) {		
				DEBUG_MSG("		Nodo " << v << " connesso a nodo ", adjacency_list[u], DEBUG_F_KERNEL);	
				DEBUG_MSG("		is_eliminated[" << adjacency_list[u] << "] -> ", is_eliminated[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		is_visited[" << adjacency_list[u] << "] -> ", is_visited[adjacency_list[u]], DEBUG_F_KERNEL);
				DEBUG_MSG("		pivots["<<v<<"] == pivots["<<adjacency_list[u]<<"] -> " << pivots[v] << " == ", pivots[adjacency_list[u]], DEBUG_F_KERNEL);

                // Si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!is_eliminated[adjacency_list[u]] && !is_visited[adjacency_list[u]] && pivots[v] == pivots[adjacency_list[u]]) {
					DEBUG_MSG("			is_visited[" << adjacency_list[u] << "] -> ", "TRUE", DEBUG_F_KERNEL);
                    // Setta il nodo puntato a visitato
					is_visited[adjacency_list[u]] = true;
                    // Permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					stop = false;
				}
			}
		}
	}
}

void reach(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_visited, bool * is_eliminated, bool * is_expanded) {
	// Esecuzione ricorsiva della chiusura in avanti
	// @param:	pivots			=	Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming

    // Tutti i pivot vengono segnati come visitati
    for(int i=0; i < num_nodes; i++) {
        is_visited[ pivots[i] ] = true;
    }

    // si effettua la chiusura in avanti
    bool stop = false;
    while(!stop) {
        stop = true;
        reach_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, is_visited, is_eliminated, is_expanded, stop);
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG("is_visited[" << i << "] -> ", is_visited[i], DEBUG_REACH);
		}
    }
}

void update(int num_nodes, int * pivots, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated, bool & stop) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che, dovrebbe contenere, per ogni 'v' dice il valore del pivot della SCC, aggiornata dopo l'esecuzione di update

    int * write_id_for_pivots = (int*) malloc(4 * num_nodes * sizeof(int));
	for (int i = 0; i < 4 * num_nodes; i++){
		write_id_for_pivots[i] = -1;
	}

	int * colors = (int*) malloc(num_nodes * sizeof(int));

	// Dai paper:
	// These subgraphs are 
	// 		1) the strongly connected component with the pivot;
	// 		2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	// 		3) the subgraph given by vertices in the backward closure but not in the forward closure;
	// 		4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	// The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	// they are recursively processed in parallel with the same algorithm

	stop = true;
	int new_color;
	for(int v = 0; v < num_nodes; v++) {
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
				//DEBUG_MSG(v, " -> non eliminato, ma non visitato da fw e bw", DEBUG_UPDATE);
			}
		}

		write_id_for_pivots[colors[v]] = v;
	}
	
	for (int i = 0; i < 4 * num_nodes; i++) {
        DEBUG_MSG("write_id_for_pivots[" + to_string(i) + "] = ", write_id_for_pivots[i], DEBUG_UPDATE);
	}

	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	for (int i = 0; i < num_nodes; i++) {
		if(is_eliminated[i]){
			pivots[i] = i;
		}else{
			pivots[i] = write_id_for_pivots[colors[i]];
		}
	}

	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_UPDATE);
	}
}

void fw_bw(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose, int *& pivots, bool * is_u) {
	// Calcola il Forward-Backward di un grafo
	// @param:	nodes						=	Lista che per ogni 'v' contiene:
	// 											nodes[v] = la poszione in 'adjacency_list' per leggere il primo nodo verso il quale parte un arco da v
	// 											nodes[v + 1] = la poszione in 'adjacency_list' per leggere l'ultimo nodo verso il quale parte un arco da v
	// 			adjacency_list				=	Lista, con cardinalità asintotica |num_edges|, che contiene il nodo ricevente di ogni arco.
	// 											È ordinata per:
	// 												1) Il numero del nodo da cui parte l'arco
	// 												2) Il numero del nodo in cui arriva l'arco (contenuto nella lista)
	// 			nodes_transpose 			=	Lista che per ogni 'u' contiene:
	// 											nodes[u] = la poszione in 'adjacency_list' per leggere il primo nodo verso il quale arriva un arco in u
	// 											nodes[u + 1] = la poszione in 'adjacency_list' per leggere l'ultimo nodo verso il quale arriva un arco in u
	// 			adjacency_list_transpose	=	Lista, con cardinalità asintotica |num_edges|, che contiene il nodo di partenza di ogni arco.
	// 											È ordinata per:
	// 												1) Il numero del nodo in cui arriva l'arco
	// 												2) Il numero del nodo da cui parte l'arco (contenuto nella lista)
	// 			pivots						= 	Lista vuota
	// 			is_u						= 	Lista che per ogni 'v' dice se fa parte di U o no
	// @return: pivots						=	Lista che per ogni 'v' dice il valore del pivot della SCC. (Le SCC possono contenere 1 solo nodo)

	bool * fw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
    bool * bw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
    bool * is_eliminated = (bool*) malloc(num_nodes * sizeof(bool));
    bool * fw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));
    bool * bw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));
    pivots = (int*) malloc(num_nodes * sizeof(int));;

	for (int i = 0; i < num_nodes; i++){
		fw_is_visited[i] = false;
		bw_is_visited[i] = false;
		is_eliminated[i] = !is_u[i];
		fw_is_expanded[i] = false;
		bw_is_expanded[i] = false;
	}

	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, is_eliminated);

	// Si prende come primo pivot globale, il primo nodo che si riesce a trovare non eliminato 
	int v = 0;
	while(v < num_nodes && is_eliminated[v]) {
		++v;
	}
	for (int i = 0; i < num_nodes; i++){
		pivots[i] = v;
	}

	DEBUG_MSG("pivot = " , v, DEBUG_FW_BW);

    bool stop = false;

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
    while (!stop){
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, nodes, adjacency_list, pivots, fw_is_visited, is_eliminated, fw_is_expanded);

		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, nodes_transpose, adjacency_list_transpose, pivots, bw_is_visited, is_eliminated, bw_is_expanded);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, nodes, nodes_transpose, adjacency_list, adjacency_list_transpose, pivots, is_eliminated);

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, pivots, fw_is_visited, bw_is_visited, is_eliminated, stop);
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

void trim_u_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_u, int *& is_scc) {
	// Setta i pivot delle SCC uguale a -1 se questi ricevono archi da nodi u
	// param: 	pivots = 	Lista che per ogni 'v' dice il valore del pivot della SCC
	// 			is_scc =	Lista copia di pivots
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1

	for(int u = 0; u < num_nodes; ++u ) {
		if(is_u[u] == true) {
			for(int v = nodes[u]; v < nodes[u+1]; ++v) {
				if(pivots[u] != pivots[adjacency_list[v]]) {
					is_scc[pivots[adjacency_list[v]]] = -1;
				}
			}
		}
	}
}

void trim_u_propagation(int num_nodes, int * pivots, int *& is_scc) {
	// Se alcuni pivot sono settati a -1, per la cancellazione dovuta a collegamenti con nodi u, 
	// propaga la cancellazione agli altri membri della SCC
	// param: 	pivots = 	Lista contenente i pivot delle SCC
	// 			is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot e gli altri nodi delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1

	for(int u = 0; u < num_nodes; ++u ) {
		is_scc[u] = is_scc[pivots[u]];
	}
}

void calculate_more_than_one(int num_nodes, int * is_scc, int *& more_than_one) {
	// Trova il numero di elementi nella SCC
	// @param: is_scc =	Lista contenente i pivot delle SCC, però i pivot e gli altri nodi delle SCC 
	// 					che ricevono archi da nodi u sono settati a -1
	// @return:	more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 								Se 'v' è pivot: 	more_than_one[v] = numero di elementi nella sua SCC,
	// 								Se 'v' non è pivot:	more_than_one[v] = 1

	for(int u = 0; u < num_nodes; ++u ) {
		if(is_scc[u] != -1)
		++more_than_one[is_scc[u]];
	}
}

void is_scc_adjust(int num_nodes, int * more_than_one, int *& is_scc) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// @param: more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 							Se 'v' è pivot: 								more_than_one[v] = numero di elementi nella sua SCC,
	// 							Se 'v' non è pivot, ma fa parte di una SCC:		more_than_one[v] = 0
	// 							Se 'v' non è pivot e non fa parte di una SCC:	more_than_one[v] = 0
	// @return: is_scc =	Lista che per ogni nodo 'v' dice se questo fa parte di una SCC.
	// 						Se fa parte di una SCC: 	is_scc[v] = valore del pivot,
	// 						Se non fa parte di una SCC:	is_scc[v] = -1

	for(int u = 0; u < num_nodes; ++u ) {
		if(more_than_one[u] == 1) {
			is_scc[u] = -1;
		}
	}
}

void trim_u(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int *& pivots, bool * is_u, int *& is_scc) {
	// Elimina le SCC riceventi archi da altri nodi U non facenti parte della SCC
	// @param:	pivots 	=	Lista che per ogni 'v' dice il valore del pivot della SCC
	// 			is_scc 	=	Lista vuota
	// 			is_u	=	Lista che per ogni 'v' dice se fa parte di U o no
	// @return: is_scc 	=	Lista che per ogni nodo 'v' dice se questo fa parte di una SCC.
	// 						Se fa parte di una SCC: 	is_scc[v] = valore del pivot,
	// 						Se non fa parte di una SCC:	is_scc[v] = -1

	is_scc = (int*) malloc(num_nodes * sizeof(int));
	for (int u = 0; u < num_nodes; ++u) {
		is_scc[u] = pivots[u];
	}

	trim_u_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, is_u, is_scc);

	trim_u_propagation(num_nodes, pivots, is_scc);

	int * more_than_one = (int*) calloc(num_nodes, sizeof(int));
	memset(more_than_one, 0, num_nodes);
	calculate_more_than_one(num_nodes, is_scc, more_than_one);

	is_scc_adjust(num_nodes, more_than_one, is_scc);
}

int count_distinct(int arr[], int n){
	// Conta quanti elementi distinti ci sono in un array
	// @param:	arr =	Array in cui contare il numero di elementi diverso
	// 			n 	=	Numero di elementi nell'array
	// @return:	res =	Numero di elementi diversi nell'array

    int res = 0;
 
    // Pick all elements one by one
    for (int i = 1; i < n; i++) {
        int j = 0;
        for (j = 0; j < i; j++)
            if (arr[i] == arr[j])
                break;
 
        // If not printed earlier, then print it
        if (i == j)
            res++;
    }
    return res;
}

int main(int argc, char ** argv) {
    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	int num_nodes, num_edges;
    int * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose, * pivots, * is_scc;
	bool * is_u;

    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, is_u);

	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_MAIN);
	}
	for (int i = 0; i < num_edges; i++) {
        DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_MAIN);
	}
	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_MAIN);
	}
	for (int i = 0; i < num_edges; i++) {
        DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_MAIN);
	}
	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("is_u[" + to_string(i) + "] = ", is_u[i], DEBUG_MAIN);
	}

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, pivots, is_u);

	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_MAIN);
	}

	trim_u(num_nodes, num_edges, nodes, adjacency_list, pivots, is_u, is_scc);

	for (int i = 0; i < num_nodes; i++) {
        DEBUG_MSG("is_scc[" + to_string(i) + "] = ", is_scc[i], false);
	}

	DEBUG_MSG("Number of SCCs found: ", count_distinct(is_scc, num_nodes), true);
}