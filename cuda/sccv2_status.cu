#include "../utils/is_checked.cu"
#include "../utils/file2graph.cpp"
#include "scc_operations.cu"
#include <cstring>
#include <cuda.h>
#include <set>
#include <unordered_map>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE false
#define DEBUG_FW_BW false
#define DEBUG_MAIN false

#ifndef DEBUG_FINAL
#define DEBUG_FINAL true
#endif

void reach_v2(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *), bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Esecuzione ricorsiva della chiusura in avanti/indietro
	// @param:	pivots			=	Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming

	*stop = false;

    // Si effettua la chiusura in avanti/indietro
    while(!*stop) {
		HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
        f_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_stop, get_visited, get_expanded, set_visited, set_expanded);
		HANDLE_ERROR(cudaMemcpy(stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }

	*stop = false;
	HANDLE_ERROR(cudaMemset(d_stop, false, sizeof(bool)));
}

void trimming_v2(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list, unsigned int * d_adjacency_list_transpose, char * d_status, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Elimina iterativamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming

	*stop = false;
    while(!*stop) {
		HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
        trimming_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, d_stop);
		HANDLE_ERROR(cudaMemcpy(stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }

	*stop = false;
	HANDLE_ERROR(cudaMemset(d_stop, false, sizeof(bool)));
}

void update_v2(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status,  unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, aggiornata dopo l'esecuzione di update
	
	*stop = true;
	HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
	// Dai paper:
	// These subgraphs are 
	// 		1) the strongly connected component with the pivot;
	// 		2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	// 		3) the subgraph given by vertices in the backward closure but not in the forward closure;
	// 		4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	// The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	// they are recursively processed in parallel with the same algorithm
	
	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	set_colors<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots, d_stop);
	cudaDeviceSynchronize();
	
	HANDLE_ERROR(cudaMemcpy(stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemset(d_stop, false, sizeof(bool)));
	
	set_new_pivots<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
	cudaDeviceSynchronize();
	set_new_eliminated<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
	cudaDeviceSynchronize();
}

void routine_v2(unsigned int num_nodes, unsigned int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

	bool * d_stop;
	bool stop;
	stop = false;

	// Dichiarazioni di variabili device
	unsigned int * d_nodes, * d_adjacency_list, * d_nodes_transpose, * d_adjacency_list_transpose, * d_pivots, * d_colors;
	char * d_status;
	unsigned long * d_write_id_for_pivots;

	const unsigned int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const unsigned int NUMBER_OF_BLOCKS = num_nodes / THREADS_PER_BLOCK + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);
	const unsigned int NUMBER_OF_BLOCKS_VEC_ACC = min(((num_nodes/4 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), prop.maxGridSize[1]);

	// Inizializzazione e copia delle funzioni device che verranno passate tramite parametro.
	// Utilizzando le funzioni in questo modo, anche se apparentemente verboso, permette di ottenere meno codice duplicato:
	// infatti, se non fosse per queste variabili, si sarebbe dovuto duplicare l'f_kernel e il reach per averne uno per la forward e uno per la backward.
	get_status h_get_fw_visited, h_get_bw_visited, h_get_fw_expanded, h_get_bw_expanded;
	set_status h_set_fw_visited, h_set_bw_visited, h_set_fw_expanded, h_set_bw_expanded;

	HANDLE_ERROR(cudaMalloc((void**)&d_write_id_for_pivots, 4 * num_nodes * sizeof(unsigned long)));
	HANDLE_ERROR(cudaMalloc((void**)&d_colors, num_nodes * sizeof(unsigned int)));

	HANDLE_ERROR(cudaMalloc((void**)&d_adjacency_list, num_edges * sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_adjacency_list_transpose, num_edges * sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_nodes, (num_nodes+1) * sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_nodes_transpose, (num_nodes+1) * sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pivots, num_nodes * sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_status, (num_nodes+1) * sizeof(char)));
	HANDLE_ERROR(cudaMalloc((void**)&d_stop, sizeof(bool)));

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpy(d_adjacency_list, adjacency_list, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_nodes, nodes, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_status, status, (num_nodes+1) * sizeof(char), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemset(d_stop, false, sizeof(bool)));

	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_fw_visited, dev_get_fw_visited, sizeof(get_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_bw_visited, dev_get_bw_visited, sizeof(get_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_fw_expanded, dev_get_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_bw_expanded, dev_get_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault));
	
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_fw_visited, dev_set_fw_visited, sizeof(set_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_bw_visited, dev_set_bw_visited, sizeof(set_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_fw_expanded, dev_set_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_bw_expanded, dev_set_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault));
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming_v2(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
	
	// Sincronizzazione implicita perché si utilizza il default stream
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();
	set_initialize_pivot<<<NUMBER_OF_BLOCKS_VEC_ACC, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();

	unsigned int * pivot_test = (unsigned int*) malloc(num_nodes * sizeof(unsigned int));
	char * status_test = (char*) malloc(num_nodes);

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
	stop = false;
    while (!stop){
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach_v2(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach_v2(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_status, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update_v2(num_nodes, d_pivots, d_status, d_colors, d_write_id_for_pivots, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		cudaDeviceSynchronize();

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		if(!stop){
			DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
			trimming_v2(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		}
    }

	HANDLE_ERROR(cudaFree(d_write_id_for_pivots));
	HANDLE_ERROR(cudaFree(d_colors));
	
	// Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1	
	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status);

	HANDLE_ERROR(cudaFree(d_adjacency_list_transpose));
	HANDLE_ERROR(cudaFree(d_adjacency_list));
	HANDLE_ERROR(cudaFree(d_nodes_transpose));
	HANDLE_ERROR(cudaFree(d_nodes));
	
	bool * d_is_scc;
	HANDLE_ERROR(cudaMalloc((void**)&d_is_scc, num_nodes * sizeof(unsigned int)));
	trim_u_propagation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status, d_is_scc);

	eliminate_trivial_scc<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(unsigned int) + THREADS_PER_BLOCK*sizeof(bool)>>>(THREADS_PER_BLOCK, num_nodes, d_pivots, d_is_scc);
	
	// Si prende come pivot, il primo pivot che si riesce a trovare facente parte di una scc 
	unsigned pivot_riferimento;
	bool pivot_riferimento_found = false;
	unsigned int * d_pivots_riferimento;
	bool * d_pivots_riferimento_found;

	HANDLE_ERROR(cudaMalloc((void**)&d_pivots_riferimento, sizeof(unsigned int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pivots_riferimento_found, sizeof(bool)));
	HANDLE_ERROR(cudaMemcpy(d_pivots_riferimento_found, &pivot_riferimento_found, sizeof(bool), cudaMemcpyHostToDevice));

	choose_scc_to_print<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_is_scc, d_pivots, d_pivots_riferimento_found, d_pivots_riferimento);
	
	HANDLE_ERROR(cudaMemcpy(&pivot_riferimento_found, d_pivots_riferimento_found, sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&pivot_riferimento, d_pivots_riferimento, sizeof(unsigned int), cudaMemcpyDeviceToHost));

	if (pivot_riferimento_found){
		if (DEBUG_FINAL){
			print_scc<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, pivot_riferimento);
			printf("\n");
		}
	}else{
		DEBUG_MSG("No SCCs found", "", DEBUG_FINAL);
	}

	HANDLE_ERROR(cudaFree(d_pivots_riferimento_found));
	HANDLE_ERROR(cudaFree(d_pivots_riferimento));
	HANDLE_ERROR(cudaFree(d_is_scc));
	HANDLE_ERROR(cudaFree(d_status));
	HANDLE_ERROR(cudaFree(d_pivots));
	HANDLE_ERROR(cudaFree(d_stop));
	
	/* cudaFreeHost(h_get_fw_visited);
	cudaFreeHost(h_get_bw_visited);
	cudaFreeHost(h_set_fw_visited);
	cudaFreeHost(h_set_bw_visited);

	cudaFreeHost(h_get_fw_expanded);
	cudaFreeHost(h_get_bw_expanded);
	cudaFreeHost(h_set_fw_expanded);
	cudaFreeHost(h_set_bw_expanded); */
}