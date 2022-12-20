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
#define DEBUG_FINAL true

#define CUDA_STREAMS 9

/*

VERSIONE DEL CODICE CUDA: SCCv4 - Pinned Memory

Questa versione del codice è un miglioramento della versione naive, in quanto si è andato a ottimizzare molti aspetti del codice:
- Le operazioni sulla memoria adesso vengono eseguite su stream diversi, sincronizzando il codice quando necessario
- Creazione di un'unica variabile "stop" da usare nei vari passaggi principali: si evita ogni volta una nuova allocazione
- Utilizzo dei registri all'interno dei kernel, per velocizzare le operazioni
- Utilizzo di un doppio shift, rimipazzando la moltiplicazione per 4
- set_colors e set_race_winners sono stati uniti, evitando il lancio di un kernel non essenziale
- Rimozione dell'array colors: tramite la programmazione parallela possiamo usare una sola variabile
- Alcune operazioni binarie su "status" sono state unite in una sola (es. 100 | 010 | 001 === 110 | 001 )
- Utilizzo di unsigned int. Nonostante non vengano migliorate direttamente le performace, c'è la possibilità di poter elaborare un numero più alto di nodi/archi
  senza dove aumentare lo spazio utilzzato in memoria

N.B.
- Non è possibile fare uso della memoria shared visto il tipo operazioni eseguite: spesso, ad esempio tramite i pivots, il codice "salta"
  da una posizione all'altra per accedere ai vari nodi. Visto che il massimo possibile sarebbe solo di salvare un frammento delle liste in memoria shared ed
  è possibile anticipare dove si andrà a leggere gli array, è impossibile farne uso.
- Non è possibile l'esecuzione in contemporanea (su stream diversi) di kernel diversi. Per funzionare correttamente ogni kernel deve ricevere i risultati di quello prima.
  L'unico caso che sarebbe possibile è quello del forward e backward reach, se non fosse che entrambi modificano l'array "status" e ci sarebbe una race condition non favorevole

*/

void reach_v4(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *), bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
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
		*stop = true;
        f_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_stop, get_visited, get_expanded, set_visited, set_expanded);
		cudaDeviceSynchronize();
    }
}

void trimming_v4(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list, unsigned int * d_adjacency_list_transpose, char * d_status, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
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

void update_v4(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status,  unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * stop, bool * d_stop, const unsigned int n_blocks, const unsigned int t_per_blocks) {
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
	
	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	set_colors<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots, d_stop);
	cudaDeviceSynchronize();
	set_new_pivots<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
	cudaDeviceSynchronize();
	/* HANDLE_ERROR(cudaMemcpy(main_stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemset(d_stop, false, sizeof(bool))); */
}

void routine_v4(const bool profiling, unsigned int num_nodes, unsigned int num_edges, unsigned * nodes, unsigned * adjacency_list, unsigned * nodes_transpose, unsigned * adjacency_list_transpose, char * status) {
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	cudaSetDeviceFlags(cudaDeviceMapHost);

	bool * d_stop, * stop;

	// Dichiarazioni di variabili device
	unsigned int * d_nodes, * d_adjacency_list, * d_nodes_transpose, * d_adjacency_list_transpose, * d_pivots, * d_colors;
	char * d_status;
	unsigned long * d_write_id_for_pivots;

	HANDLE_ERROR(cudaHostRegister(nodes, (num_nodes+1) * sizeof(unsigned int), cudaHostRegisterDefault));
	HANDLE_ERROR(cudaHostRegister(adjacency_list, num_edges * sizeof(unsigned int), cudaHostRegisterDefault));
	HANDLE_ERROR(cudaHostRegister(nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaHostRegisterDefault));
	HANDLE_ERROR(cudaHostRegister(adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaHostRegisterDefault));
	HANDLE_ERROR(cudaHostRegister(status, (num_nodes+1) * sizeof(char), cudaHostRegisterDefault));

	/* if(DEBUG_MAIN){
		for (unsigned int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_MAIN);
		for (unsigned int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_MAIN);
		for (unsigned int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_MAIN);
		for (unsigned int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_MAIN);
	} */

	const unsigned int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const unsigned int NUMBER_OF_BLOCKS = num_nodes / THREADS_PER_BLOCK + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);

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

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list, adjacency_list, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[1]));
	HANDLE_ERROR(cudaMemcpyAsync(d_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[2]));
	HANDLE_ERROR(cudaMemcpyAsync(d_nodes, nodes, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[3]));
	HANDLE_ERROR(cudaMemcpyAsync(d_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(unsigned int), cudaMemcpyHostToDevice, stream[4]));
	HANDLE_ERROR(cudaMemcpyAsync(d_status, status, (num_nodes+1) * sizeof(char), cudaMemcpyHostToDevice, stream[6]));

	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_visited, dev_get_fw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[0]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_visited, dev_get_bw_visited, sizeof(get_status), 0, cudaMemcpyDefault, stream[5]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_fw_expanded, dev_get_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[1]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_get_bw_expanded, dev_get_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[2]));
	
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_visited, dev_set_fw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[3]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_visited, dev_set_bw_visited, sizeof(set_status), 0, cudaMemcpyDefault, stream[4]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_fw_expanded, dev_set_fw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[5]));
	HANDLE_ERROR(cudaMemcpyFromSymbolAsync(&h_set_bw_expanded, dev_set_bw_expanded, sizeof(get_status), 0, cudaMemcpyDefault, stream[6]));

	HANDLE_ERROR(cudaHostAlloc(&stop, sizeof(bool), cudaHostAllocMapped));
	HANDLE_ERROR(cudaHostGetDevicePointer(&d_stop, stop, 0));

	cudaStreamSynchronize(stream[1]);
	cudaStreamSynchronize(stream[2]);
	cudaStreamSynchronize(stream[3]);
	cudaStreamSynchronize(stream[4]);
	cudaStreamSynchronize(stream[6]);
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	*stop = false;
    while(!*stop) {
		*stop = true;
        trimming_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[1]>>>(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, d_stop);
    }
	
	// Sincronizzazione implicita perché si utilizza il default stream
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	cudaDeviceSynchronize();
	set_initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);

	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
	*stop = false;
    while (!*stop){
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach_v4(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach_v4(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_status, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		//DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        //trimming_v4(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update_v4(num_nodes, d_pivots, d_status, d_colors, d_write_id_for_pivots, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		if(!*stop){
			DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
			trimming_v4(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, stop, d_stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
			*stop = false;
		}
    }

	cudaFreeHost(stop);
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

	if(profiling){
		eliminate_trivial_scc<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK*sizeof(unsigned int) + THREADS_PER_BLOCK*sizeof(bool)>>>(THREADS_PER_BLOCK, num_nodes, d_pivots, d_is_scc);
		cudaDeviceSynchronize();
		
		bool result = is_there_an_scc(NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, num_nodes, d_is_scc);
		printf("%d\n", result);
	}else{
		// Nella versione naive, una funzione calcolava il numero di nodi di una SCC e poi "cancellava" quelli con un numero < 2.
		// La funzione è stata eliminata e is_scc_adjust si occupa di "cancellare" tali nodi senza doverli contare.
		// N.B. Per "cancellare" si intende assegnare ad un generico nodo v is_scc[v] = -1
		is_scc_adjust<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, 0, stream[6]>>>(num_nodes, d_pivots, d_status);
		cudaDeviceSynchronize();

		// Questa sezione di codice è temporanea, verrà rimossa al momento del test
		unsigned int * pivots;
		char * final_status;

		pivots = (unsigned int*) malloc(num_nodes * sizeof(unsigned int));
		final_status = (char*) malloc(num_nodes * sizeof(char));

		cudaMemcpy(pivots, d_pivots, num_nodes * sizeof(unsigned int), cudaMemcpyDeviceToHost);
		cudaMemcpy(final_status, d_status, num_nodes * sizeof(char), cudaMemcpyDeviceToHost);

		DEBUG_MSG("Number of SCCs found: ", count_distinct_scc(num_nodes, pivots, final_status), DEBUG_FINAL);

		HANDLE_ERROR(cudaFree(d_pivots));
		HANDLE_ERROR(cudaFree(d_status));
		free(final_status);
		free(pivots);
	}


	// Da scommentare una volta finito il progetto
	//HANDLE_ERROR(cudaFreeAsync(d_pivots, stream[0]));
	//HANDLE_ERROR(cudaFreeAsync(d_status, stream[1]));

	cudaFreeHost(h_get_fw_visited);
	cudaFreeHost(h_get_bw_visited);
	cudaFreeHost(h_set_fw_visited);
	cudaFreeHost(h_set_bw_visited);

	cudaFreeHost(h_get_fw_expanded);
	cudaFreeHost(h_get_bw_expanded);
	cudaFreeHost(h_set_fw_expanded);
	cudaFreeHost(h_set_bw_expanded);

	HANDLE_ERROR(cudaHostUnregister(nodes));
	HANDLE_ERROR(cudaHostUnregister(nodes_transpose));
	HANDLE_ERROR(cudaHostUnregister(adjacency_list));
	HANDLE_ERROR(cudaHostUnregister(adjacency_list_transpose));
	HANDLE_ERROR(cudaHostUnregister(status));

	cudaDeviceSynchronize();

	for (unsigned int i=0; i<CUDA_STREAMS; ++i){
		cudaStreamDestroy(stream[i]);
	}
}