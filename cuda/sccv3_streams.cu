#include "../utils/is_checked.cu"
#include "../utils/file2graph.cpp"
#include "scc_operations.cu"
#include <cstring>
#include <cuda.h>
#include <set>
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

#define CUDA_STREAMS 9

void reach_v3(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *), bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
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

void trimming_v3(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list, unsigned int * d_adjacency_list_transpose, char * d_status, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
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

void update_v3(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status,  unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, aggiornata dopo l'esecuzione di update
	
	*stop = true;
	HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
	
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

void routine_v3(unsigned int num_nodes, unsigned int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
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

	cudaStream_t stream[CUDA_STREAMS];
	for (unsigned int i=0; i<CUDA_STREAMS; ++i){
		cudaStreamCreate(&stream[i]);
	}

	HANDLE_ERROR(cudaMallocAsync((void**)&d_write_id_for_pivots, 4 * num_nodes * sizeof(unsigned long), stream[0]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_colors, num_nodes * sizeof(unsigned int), stream[0]));

	HANDLE_ERROR(cudaMallocAsync((void**)&d_adjacency_list, num_edges * sizeof(unsigned int), stream[1]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_adjacency_list_transpose, num_edges * sizeof(unsigned int), stream[2]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_nodes, (num_nodes+1) * sizeof(unsigned int), stream[3]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_nodes_transpose, (num_nodes+1) * sizeof(unsigned int), stream[4]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_pivots, num_nodes * sizeof(unsigned int), stream[5]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_status, (num_nodes+1) * sizeof(char), stream[6]));
	HANDLE_ERROR(cudaMallocAsync((void**)&d_stop, sizeof(bool), stream[7]));

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list, adjacency_list, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[1]));
	HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[2]));
	HANDLE_ERROR(cudaMemcpyAsync(d_nodes, nodes, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[3]));
	HANDLE_ERROR(cudaMemcpyAsync(d_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[4]));
	HANDLE_ERROR(cudaMemcpyAsync(d_status, status, (num_nodes+1) * sizeof(char), cudaMemcpyHostToDevice, stream[6]));
	HANDLE_ERROR(cudaMemcpyAsync(d_stop, &stop, sizeof(bool), cudaMemcpyHostToDevice, stream[7]));

	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_visited, dev_get_fw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[0]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_visited, dev_get_bw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[5]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_expanded, dev_get_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[1]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_expanded, dev_get_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[2]));
	
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_visited, dev_set_fw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[3]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_visited, dev_set_bw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[4]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_expanded, dev_set_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[5]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_expanded, dev_set_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[6]));

	cudaStreamSynchronize(stream[1]);
	cudaStreamSynchronize(stream[2]);
	cudaStreamSynchronize(stream[3]);
	cudaStreamSynchronize(stream[4]); 
	cudaStreamSynchronize(stream[6]);
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming_v3(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
	
	// Sincronizzazione implicita perché si utilizza il default stream
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();
	set_initialize_pivot<<<NUMBER_OF_BLOCKS_VEC_ACC, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
	stop = false;
    while (!stop){
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach_v3(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach_v3(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_status, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update_v3(num_nodes, d_pivots, d_status, d_colors, d_write_id_for_pivots, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		if(!stop){
			DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
			trimming_v3(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, &stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		}
    }

	cudaFreeHost(d_stop);
	HANDLE_ERROR(cudaFreeAsync(d_write_id_for_pivots, stream[0]));
	HANDLE_ERROR(cudaFreeAsync(d_colors, stream[0]));
	
	// Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1
	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[6]>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status);
	
	HANDLE_ERROR(cudaFreeAsync(d_adjacency_list_transpose, stream[1]));
	HANDLE_ERROR(cudaFreeAsync(d_adjacency_list, stream[2]));
	HANDLE_ERROR(cudaFreeAsync(d_nodes_transpose, stream[3]));
	HANDLE_ERROR(cudaFreeAsync(d_nodes, stream[4]));
	
	cudaStreamSynchronize(stream[6]);

	bool * d_is_scc;
	HANDLE_ERROR(cudaMalloc((void**)&d_is_scc, num_nodes * sizeof(unsigned int)));
	trim_u_propagation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[6]>>>(num_nodes, d_pivots, d_status, d_is_scc);

	cudaStreamSynchronize(stream[6]);

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

	// Da scommentare una volta finito il progetto
	HANDLE_ERROR(cudaFreeAsync(d_pivots, stream[0]));
	HANDLE_ERROR(cudaFreeAsync(d_status, stream[1]));
	HANDLE_ERROR(cudaFreeAsync(d_is_scc, stream[2]));

	HANDLE_ERROR(cudaFree(d_pivots_riferimento_found));
	HANDLE_ERROR(cudaFree(d_pivots_riferimento));
/* 	cudaFreeHost(h_get_fw_visited);
	cudaFreeHost(h_get_bw_visited);
	cudaFreeHost(h_set_fw_visited);
	cudaFreeHost(h_set_bw_visited);

	cudaFreeHost(h_get_fw_expanded);
	cudaFreeHost(h_get_bw_expanded);
	cudaFreeHost(h_set_fw_expanded);
	cudaFreeHost(h_set_bw_expanded); */

	cudaDeviceSynchronize();

	for (unsigned int i=0; i<CUDA_STREAMS; ++i){
		cudaStreamDestroy(stream[i]);
	}
}