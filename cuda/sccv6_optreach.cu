#include "../utils/is_checked.cu"
#include "../utils/file2graph.cpp"
#include "scc_operations.cu"
#include <cstring>
#include <cuda.h>
#include <omp.h>
#include <set>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE false
#define DEBUG_FW_BW false
#define DEBUG_MAIN false
#define DEBUG_FINAL true

#define CUDA_STREAMS 9
#define OMP_MIN_NODES 100000
/*
VERSIONE DEL CODICE CUDA: SCCv5 - OpenMP
Rispetto alla quarta versione, in questa vengono parallelizzati, tramite le direttive di OpenMP, le chiamate all'API di CUDA per l'allocazione e il trasferimento dei dati
*/

void trimming_v6(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list, unsigned int * d_adjacency_list_transpose, char * d_status, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Elimina iterativamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming

	*stop = false;
    while(!*stop) {
		*stop = true;
        trimming_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, d_stop);
		// Dobbiamo aspettare che la copia venga effettuata anche se è mappata
		cudaDeviceSynchronize();
    }
}

void update_v6(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status, unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, aggiornata dopo l'esecuzione di update
	
	*d_stop = true;

	// Dai paper:
	// These subgraphs are 
	// 		1) the strongly connected component with the pivot;
	// 		2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	// 		3) the subgraph given by vertices in the backward closure but not in the forward closure;
	// 		4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	// The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	// they are recursively processed in parallel with the same algorithm
	
	set_colors<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots, d_stop);
	cudaDeviceSynchronize();

	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	set_new_pivots<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
}

void routine_v6(const bool profiling, unsigned int num_nodes, unsigned int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	cudaSetDeviceFlags(cudaDeviceMapHost);

	const short MAX_THREADS_OMP = omp_get_max_threads();

	bool * stop, * d_stop, * bw_stop, * d_bw_stop;

	// Dichiarazioni di variabili device
	unsigned int * d_nodes, * d_adjacency_list, * d_nodes_transpose, * d_adjacency_list_transpose, * d_pivots, * d_colors;
	char * d_status, * d_bw_status;
	unsigned long * d_write_id_for_pivots;

	// Page-locking delle strutture dati principali
	#pragma omp parallel sections if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{
		#pragma omp section 
		{
			HANDLE_ERROR(cudaHostRegister(nodes, (num_nodes+1) * sizeof(unsigned int), cudaHostRegisterDefault));
		}

		#pragma omp section 
		{
			HANDLE_ERROR(cudaHostRegister(adjacency_list, num_edges * sizeof(unsigned int), cudaHostRegisterDefault));
		}

		#pragma omp section 
		{
			HANDLE_ERROR(cudaHostRegister(nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaHostRegisterDefault));
		}

		#pragma omp section 
		{
			HANDLE_ERROR(cudaHostRegister(adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaHostRegisterDefault));
		}

		#pragma omp section 
		{
			HANDLE_ERROR(cudaHostRegister(status, (num_nodes+1) * sizeof(char), cudaHostRegisterDefault));
		}
	}

	const unsigned int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const unsigned int NUMBER_OF_BLOCKS = (num_nodes / THREADS_PER_BLOCK) + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);

	// Inizializzazione e copia delle funzioni device che verranno passate tramite parametro.
	// Utilizzando le funzioni in questo modo, anche se apparentemente verboso, permette di ottenere meno codice duplicato:
	// infatti, se non fosse per queste variabili, si sarebbe dovuto duplicare l'f_kernel e il reach per averne uno per la forward e uno per la backward.
	get_status h_get_fw_visited, h_get_bw_visited, h_get_fw_expanded, h_get_bw_expanded;
	set_status h_set_fw_visited, h_set_bw_visited, h_set_fw_expanded, h_set_bw_expanded;

	cudaStream_t stream[CUDA_STREAMS];

	// Creazione delle stream, allocazione delle variabili device e copia dei dati
	// La parallelizzazione della copia dei dati non è effettiva, in quanto il canale di comunicazione è solo uno
	// e quindi la copia dei dati è serializzata.
	#pragma omp parallel if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{
		#pragma omp for schedule(static) 
		for (short i = 0; i < CUDA_STREAMS; i++) {
			cudaStreamCreate(&stream[i]);
		}

		#pragma omp barrier

		#pragma omp sections nowait
		{
			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_visited, dev_get_fw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[0]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_visited, dev_get_bw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[1]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_expanded, dev_set_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[2]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_expanded, dev_get_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[3]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_expanded, dev_get_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[4]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_visited, dev_set_fw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[5]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_visited, dev_set_bw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[6]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_expanded, dev_set_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[7]));
			}	
		}

		//Allocazione delle variabili device
		#pragma omp sections nowait
		{
			#pragma omp section 
			{
				HANDLE_ERROR(cudaMallocAsync((void**)&d_write_id_for_pivots, 4 * num_nodes * sizeof(unsigned long), stream[0]));
				HANDLE_ERROR(cudaMallocAsync((void**)&d_pivots, num_nodes * sizeof(unsigned int), stream[1]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMallocAsync((void**)&d_adjacency_list, num_edges * sizeof(unsigned int), stream[2]));
				HANDLE_ERROR(cudaMallocAsync((void**)&d_adjacency_list_transpose, num_edges * sizeof(unsigned int), stream[3]));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaMallocAsync((void**)&d_nodes, (num_nodes+1) * sizeof(unsigned int), stream[4]));
				HANDLE_ERROR(cudaMallocAsync((void**)&d_nodes_transpose, (num_nodes+1) * sizeof(unsigned int), stream[5]));
			}
			
			#pragma omp section 
			{
				HANDLE_ERROR(cudaMallocAsync((void**)&d_status, (num_nodes+1) * sizeof(char), stream[6]));
				HANDLE_ERROR(cudaMallocAsync((void**)&d_bw_status, (num_nodes+1) * sizeof(char), stream[7]));
				HANDLE_ERROR(cudaMallocAsync((void**)&d_colors, num_nodes * sizeof(unsigned int), stream[8]));
			}
		}

		#pragma omp barrier
		
		// Sincronizzazione delle stream
		#pragma omp for schedule(static)
		for (short i = 2; i < 7; i++) {
			cudaStreamSynchronize(stream[i]);
		}

		#pragma omp barrier

		// cudaMemcpy per archi e nodi
		#pragma omp sections
		{
			#pragma omp section
			{
				HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list, adjacency_list, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[2]));				
			}

			#pragma omp section 
			{	
				HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[3]));
			}

			#pragma omp section
			{
				HANDLE_ERROR(cudaMemcpyAsync(d_nodes, nodes, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[4]));				
			}

			#pragma omp section
			{
				HANDLE_ERROR(cudaMemcpyAsync(d_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[5]));	
			} 

			#pragma omp section
			{
				HANDLE_ERROR(cudaMemcpyAsync(d_status, status, (num_nodes+1) * sizeof(char), cudaMemcpyHostToDevice, stream[6]));				
			}
		}	
	}

	HANDLE_ERROR(cudaHostAlloc(&stop, sizeof(bool), cudaHostAllocMapped));
	HANDLE_ERROR(cudaHostGetDevicePointer(&d_stop, stop, 0));

	HANDLE_ERROR(cudaHostAlloc(&bw_stop, sizeof(bool), cudaHostAllocMapped));
	HANDLE_ERROR(cudaHostGetDevicePointer(&d_bw_stop, bw_stop, 0));
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming_v6(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

	/* Inizializzazione delle variabili di test per il controllo della correttezza dell'algoritmo
	char * status_tmp, * bw_status_tmp;
	unsigned int * pivots_tmp;
	pivots_tmp = (unsigned int *) malloc(num_nodes * sizeof(unsigned int));
	status_tmp = (char *) malloc(num_nodes * sizeof(char));
	bw_status_tmp = (char *) malloc(num_nodes * sizeof(char)); */

	// Sincronizzazione implicita perché si utilizza il default stream
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();
	set_initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	
	/* Print di debug riguardante lo stato dei nodi e i pivot iniziali
	
	HANDLE_ERROR(cudaMemcpy(status_tmp, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(pivots_tmp, d_pivots, num_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost));
	for(int i = 0; i < num_nodes; i++){
		printf("status[%d] = %s, pivots[%d] = %d\n", i, from_status_to_string(status_tmp[i]), i, pivots_tmp[i]);
	} */

	HANDLE_ERROR(cudaMemcpy(d_bw_status, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToDevice));

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
	*stop = false;
	*bw_stop = false;
    while (!*stop){
		// Forward + Backward reach
		DEBUG_MSG("Reach:" , "", DEBUG_FW_BW);
       	
		*stop = false;
		*bw_stop = false;
		HANDLE_ERROR(cudaMemcpy(d_bw_status, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToDevice));
		while(!(*stop && *bw_stop)) {
			
			#pragma omp parallel sections if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
			{
				#pragma omp section 
				{
					if(!*stop){
						*stop = true;
						f_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[1]>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_stop, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded);
					}
				}
				#pragma omp section 
				{
					if(!*bw_stop){
						*bw_stop = true;
						f_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[2]>>>(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_bw_status, d_bw_stop, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded);
					}
				}
			}

			cudaStreamSynchronize(stream[1]);
			cudaStreamSynchronize(stream[2]);
		}

		/*  Print di debug riguardante lo stato dei nodi e i pivot
		
		HANDLE_ERROR(cudaMemcpy(status_tmp, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(bw_status_tmp, d_bw_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost));
		for(int i = 0; i < num_nodes; i++){
			printf("fw_status[%d] = %s, bw_status[%d] = %s\n", i, from_status_to_string(status_tmp[i]), i, from_status_to_string(bw_status_tmp[i]));
		}
		for(int i = 0; i < num_nodes; i++){
			printf("pivots[%d] = %d\n", i, pivots_tmp[i]);
		} */
 		
		bitwise_or_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_status, d_bw_status);
		cudaDeviceSynchronize();

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update_v6(num_nodes, d_pivots, d_status, d_colors, d_write_id_for_pivots, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		if(!*stop){
			DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
			trimming_v6(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
			*stop = false;
		}

		/* Print riguardante i nodi elimitati
		HANDLE_ERROR(cudaMemcpy(status_tmp, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost));
		for(int i = 0; i < num_nodes; i++){
			printf("status[%d] = %s\n", i, from_status_to_string(status_tmp[i]));
		} */

		/* Print riguardante i nuovi pivot
		HANDLE_ERROR(cudaMemcpy(pivots_tmp, d_pivots, num_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost));
		for(int i = 0; i < num_nodes; i++){
			printf("pivots[%d] = %d\n", i, pivots_tmp[i]);
		}
		printf("---------------------\n");

		cudaDeviceSynchronize();
		HANDLE_ERROR(cudaMemcpy(status_tmp, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		*/
    }
	
	#pragma omp parallel sections if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{
		#pragma omp section
		{
			cudaFreeHost(stop);
			cudaFreeHost(d_stop);
		}

		#pragma omp section
		{
			HANDLE_ERROR(cudaFreeAsync(d_write_id_for_pivots, stream[0]));
		}

		#pragma omp section
		{
			HANDLE_ERROR(cudaFreeAsync(d_colors, stream[1]));
		}

		#pragma omp section
		{
			cudaFreeHost(h_get_fw_visited);
			cudaFreeHost(h_get_bw_visited);
		}

		#pragma omp section
		{
			cudaFreeHost(h_set_fw_visited);
			cudaFreeHost(h_set_bw_visited);
		}

		#pragma omp section
		{
			cudaFreeHost(h_get_fw_expanded);
			cudaFreeHost(h_get_bw_expanded);
		}

		#pragma omp section
		{
			cudaFreeHost(h_set_fw_expanded);
			cudaFreeHost(h_set_bw_expanded);
			cudaFreeHost(h_set_bw_expanded);
		}
	}

	// Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1
	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status);
	
	bool * d_is_scc;
    #pragma omp parallel sections if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{	
		#pragma omp section
		{
			HANDLE_ERROR(cudaHostUnregister(adjacency_list_transpose));
			HANDLE_ERROR(cudaFreeAsync(d_adjacency_list_transpose, stream[1]));
		}

		#pragma omp section
		{
			HANDLE_ERROR(cudaHostUnregister(adjacency_list));
			HANDLE_ERROR(cudaFreeAsync(d_adjacency_list, stream[2]));
		}

		#pragma omp section
		{
			HANDLE_ERROR(cudaHostUnregister(nodes_transpose));
			HANDLE_ERROR(cudaFreeAsync(d_nodes_transpose, stream[3]));
		}

		#pragma omp section
		{
			HANDLE_ERROR(cudaHostUnregister(nodes));
			HANDLE_ERROR(cudaFreeAsync(d_nodes, stream[4]));
		}

		#pragma omp section
		{	
			HANDLE_ERROR(cudaMalloc((void**)&d_is_scc, num_nodes * sizeof(unsigned int)));
			trim_u_propagation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status, d_is_scc);
		}
    }

	if(profiling){
		eliminate_trivial_scc<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(unsigned int) + THREADS_PER_BLOCK*sizeof(bool)>>>(THREADS_PER_BLOCK, num_nodes, d_pivots, d_is_scc);
		
		bool result = is_there_an_scc(NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, num_nodes, d_is_scc);
		printf("%d\n", result);
	}else{
		// Nella versione naive, una funzione calcolava il numero di nodi di una SCC e poi "cancellava" quelli con un numero < 2.
		// La funzione è stata eliminata e is_scc_adjust si occupa di "cancellare" tali nodi senza doverli contare.
		// N.B. Per "cancellare" si intende assegnare ad un generico nodo v is_scc[v] = -1
		is_scc_adjust<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
		cudaDeviceSynchronize();
		is_scc_adjust_prop<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);

		unsigned int * pivots;
		char * final_status;

		pivots = (unsigned int*) malloc(num_nodes * sizeof(unsigned int));
		final_status = (char*) malloc(num_nodes * sizeof(char));

		cudaMemcpy(pivots, d_pivots, num_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost);
		cudaMemcpy(final_status, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost);

		DEBUG_MSG("Number of SCCs found: ", count_distinct_scc(num_nodes, pivots, final_status), DEBUG_FINAL);

		free(final_status);
		free(pivots);
	}

	#pragma omp parallel if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{
		#pragma omp sections nowait
		{
			#pragma omp section 
			{
				HANDLE_ERROR(cudaFree(d_is_scc));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaFree(d_status));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaFree(d_pivots));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaHostUnregister(status));
			}
		}

		#pragma omp barrier

		#pragma omp for schedule(static)
		for (short i = 0; i < CUDA_STREAMS; i++) {
			cudaStreamDestroy(stream[i]);
		}
	}


}