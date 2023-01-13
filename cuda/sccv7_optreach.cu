#include "../utils/is_checked.cu"
#include "../utils/file2graph.cpp"
#include "scc_operations.cu"
#include <cstring>
#include <cuda.h>
#include <omp.h>
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

#define CUDA_STREAMS 9
#ifndef OMP_MIN_NODES
#define OMP_MIN_NODES 10000
#endif

void trimming_v7(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list, unsigned int * d_adjacency_list_transpose, char * d_status, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Elimina iterativamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati

	*stop = false;
    while(!*stop) {
		*stop = true;
        trimming_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, d_stop);
		// Dobbiamo aspettare che l'esecuzione termini prima di provare ad eseguire la prossima iterazione
		cudaDeviceSynchronize();
    }
}

void update_v7(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status, unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
	// Esegue l'update dei valori del pivot, facendo una race
	
	*d_stop = true;
	
	// set_colors andrà ad assegnare una nuova sottorete (colore) a tutti i nodi che non sono stati eliminati (quindi chi non è ancora in una SCC)
	set_colors<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots, d_stop);
	cudaDeviceSynchronize();

	// Setto i valori dei pivot che hanno vinto la race
	set_new_pivots<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
	cudaDeviceSynchronize();
	set_new_eliminated<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
}

void routine_v7(unsigned int num_nodes, unsigned int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

	// Si utilizza il flag cudaDeviceMapHost per poter faciliare il passaggio di dati della variabile stop tra host e device
	cudaSetDeviceFlags(cudaDeviceMapHost);

	// Dichiarazione del numero massimo di thread utilizzabili da OpenMP
	const short MAX_THREADS_OMP = omp_get_max_threads();
	
	// Dichiara le variabili per la terminazione
	bool * stop, * d_stop, * bw_stop, * d_bw_stop;

	// Dichiarazioni di variabili device
	unsigned int * d_nodes, * d_adjacency_list, * d_nodes_transpose, * d_adjacency_list_transpose, * d_pivots, * d_colors;
	char * d_status;
	unsigned long * d_write_id_for_pivots;

	// Page-locking delle strutture dati principali tramite OpenMP
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

	// Inizializzazione delle costanti utilizzate per la gestione dei blocchi e dei thread
	// NUMBER_OF_BLOCKS_VEC_ACC è il numero di blocchi utilizzati per l'accesso vettorizzato ai vettori
	const unsigned int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const unsigned int NUMBER_OF_BLOCKS = (num_nodes / THREADS_PER_BLOCK) + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);
	const unsigned int NUMBER_OF_BLOCKS_VEC_ACC = min(((num_nodes/4 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), prop.maxGridSize[1]);

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
		// Creazione degli stream
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

	// Mappatura della variabile stop (forward + restanti utilizzi)
	HANDLE_ERROR(cudaHostAlloc(&stop, sizeof(bool), cudaHostAllocMapped));
	HANDLE_ERROR(cudaHostGetDevicePointer(&d_stop, stop, 0));

	// Mappatura della variabile stop (utilizzata dal backward reach)
	HANDLE_ERROR(cudaHostAlloc(&bw_stop, sizeof(bool), cudaHostAllocMapped));
	HANDLE_ERROR(cudaHostGetDevicePointer(&d_bw_stop, bw_stop, 0));
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U non avevano più out-degree e in-degree diverso da 0
	trimming_v7(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

	// Sincronizzazione implicita perché si utilizza il default stream
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();
	set_initialize_pivot<<<NUMBER_OF_BLOCKS_VEC_ACC, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
	*stop = false;
	*bw_stop = false;
    while (!*stop){
		// Forward + Backward reach
		DEBUG_MSG("Reach:" , "", DEBUG_FW_BW);
       	
		*stop = false;
		*bw_stop = false;

		// Parallelizazione dei due reach
		#pragma omp parallel sections num_threads(2)
		{
			// Forward reach
			#pragma omp section 
			{
				while(!*stop){
					*stop = true;
					f_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[1]>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_stop, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded);
					cudaStreamSynchronize(stream[1]);
				}
			}
			// Backward reach
			#pragma omp section 
			{
				while(!*bw_stop){
					*bw_stop = true;
					f_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[2]>>>(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_status, d_bw_stop, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded);
					cudaStreamSynchronize(stream[2]);
				}
			}
		}

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update_v7(num_nodes, d_pivots, d_status, d_colors, d_write_id_for_pivots, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		if(!*stop){
			DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
			trimming_v7(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
			*stop = false;
		}
    }

	// L'algoritmo di identificazione delle SCC è concluso, si procede con una liberazione parziale della memoria allocata
	
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

		/* #pragma omp section
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
		} */
	}

	// Setta i pivot delle SCC come non facenti parte di una SCC se queste ricevono archi da nodi u
	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status);
	
	// d_is_scc verrà utilizzato per identificare se un nodo è una SCC
	// Potrebbe essere evitato l'utilizzo di un array di bool, ma per semplicità di implementazione è stato utilizzato
	// Oltre che per semplicità, ci permette di gestire in modo più veloce tutte le varie versioni dell'algoritmo
	bool * d_is_scc;

	// Liberazione di altra memoria (nodi, archi, pivot, status, bw_status) + trim_u_propagation
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

	// Dato che l'algoritmo Forward-Backward identifica anche i singoli nodi come SCC
	// Setto tutti i pivot come non facenti parte di una SCC, facendo attenzione a non rieseguire la funzione trim_u_propagation.
	// Quindi tutte le SCC da 1 nodo saranno eliminate, mentre le SCC con 2 o più nodi verranno ancora considerate tali.
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

	// Liberazione di memoria finale e distruzione streams
	#pragma omp parallel if(num_nodes>OMP_MIN_NODES) num_threads(MAX_THREADS_OMP)
	{
		#pragma omp sections
		{
			#pragma omp section 
			{
				HANDLE_ERROR(cudaFree(d_is_scc));
				HANDLE_ERROR(cudaFree(d_pivots_riferimento));
			}

			#pragma omp section 
			{
				HANDLE_ERROR(cudaFree(d_status));
				HANDLE_ERROR(cudaFree(d_pivots_riferimento_found));
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

		#pragma omp for schedule(static)
		for (short i = 0; i < CUDA_STREAMS; i++) {
			cudaStreamDestroy(stream[i]);
		}
	}


}