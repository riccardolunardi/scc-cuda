#include "../utils/is_checked.cu"
#include "../utils/file2graph.cpp"
#include <cstring>
#include <cuda.h>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH false
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING false
#define DEBUG_UPDATE false
#define DEBUG_FW_BW false
#define DEBUG_MAIN false
#define DEBUG_FINAL true

static void handle_error(cudaError_t err, const char *file, int line ) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
		exit( EXIT_FAILURE );
	}
}
#define HANDLE_ERROR( err ) (handle_error( err, __FILE__, __LINE__ ))

__global__ void f_kernel(int num_nodes, int * d_nodes, int * d_adjacency_list, int * d_pivots, char * d_status, bool * d_stop, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *)){
	// Esecuzione di un thread della chiusura in avanti/indietro
	// @param:	pivots			=	Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming
	
	int v = threadIdx.x + blockIdx.x * blockDim.x;

    // Per ogni nodo
	if(v < num_nodes) {
        // Si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!get_is_d_eliminated(&d_status[v]) && get_visited(&d_status[v]) && !get_expanded(&d_status[v])) {
            // Si segna come espanso
			set_expanded(&d_status[v]);

            // Per ogni nodo a cui punta
			for(int u = d_nodes[v]; u < d_nodes[v + 1]; u++) {	
				int dst = d_adjacency_list[u];

                // Si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!get_is_d_eliminated(&d_status[dst]) && !get_visited(&d_status[dst]) && d_pivots[v] == d_pivots[dst]) {
                    // Setta il nodo puntato a visitato
					set_visited(&d_status[dst]);
                    // Permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					*d_stop = false;
				}
			}
		}
	}
}

void reach(int num_nodes, int * d_nodes, int * d_adjacency_list, int * d_pivots, char * d_status, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *), const int t_per_blocks,  const int n_blocks) {
	// Esecuzione ricorsiva della chiusura in avanti/indietro
	// @param:	pivots			=	Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming

	bool stop, *d_stop;
	stop = false;

	HANDLE_ERROR(cudaMalloc((void**)&d_stop, sizeof(bool)));
	
    // Si effettua la chiusura in avanti/indietro
    while(!stop) {
		HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
        f_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_stop, get_visited, get_expanded, set_visited, set_expanded);	
		HANDLE_ERROR(cudaMemcpy(&stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }
	
	HANDLE_ERROR(cudaFree(d_stop));
}

__global__ void trimming_kernel(int num_nodes, int * d_nodes, int * d_nodes_transpose, int * d_adjacency_list,  int * d_adjacency_list_transpose, char * d_status, bool * d_stop){
	// Esegue un'eliminazione di nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming
	
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		if(!get_is_d_eliminated(&d_status[v])){
			// Se questo valore non verrà cambiato, allora il nodo verrà cancellato
			bool elim = true;

			bool forward = false;
			bool backward = false;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0, tra i soli nodi non eliminati, allora non va eliminato
			for(int u = d_nodes[v]; u < d_nodes[v+1]; u++){
				if(!get_is_d_eliminated(&d_status[d_adjacency_list[u]])) {
					forward = true;
				}
			}
			if(forward) {
				for(int u = d_nodes_transpose[v]; u < d_nodes_transpose[v+1]; u++){
					if(!get_is_d_eliminated(&d_status[d_adjacency_list_transpose[u]])) {
						backward = true;
					}
				}
			}
			if(backward) {
				elim = false;
			}

			if(elim){
				set_is_d_eliminated(&d_status[v]);
				*d_stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int * d_nodes, int * d_nodes_transpose, int * d_adjacency_list, int * d_adjacency_list_transpose, char * d_status, const int n_blocks, const int t_per_blocks) {
	// Elimina iterativamente i nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming

	bool stop, *d_stop;
	stop = false;

	HANDLE_ERROR(cudaMalloc((void**)&d_stop, sizeof(bool)));

    while(!stop) {
		HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
        trimming_kernel<<<n_blocks, t_per_blocks>>>(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, d_stop);
		HANDLE_ERROR(cudaMemcpy(&stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }

	HANDLE_ERROR(cudaFree(d_stop));
}

__global__ void set_colors(int num_nodes, char * d_status, int * d_pivots, int * d_colors, long * d_write_id_for_pivots, bool * d_stop){
	// Esegue l'update dei valori del pivot facendo una race, scrivendo il "colore" di una serie di pivot in array simultaneamente
	// @param:	pivots						= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated				= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: d_write_id_for_pivots		= Lista che conterrà, nelle posizione identificate dai colori appena calcolati, i nuovi pivot da assegnare
	
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		if(get_is_d_eliminated(&d_status[v])){
			d_pivots[v] = v;
		}

		const bool fw_visitated = get_is_d_fw_visited(&d_status[v]);
		const bool bw_visitated = get_is_d_bw_visited(&d_status[v]);
		
		if(fw_visitated == bw_visitated && fw_visitated == true){
			d_colors[v] = 4 * d_pivots[v];
		} else {
			if(fw_visitated != bw_visitated && fw_visitated == true){
				d_colors[v] = 4 * d_pivots[v] + 1;
			}else if(fw_visitated != bw_visitated && fw_visitated == false){
				d_colors[v] = 4 * d_pivots[v] + 2;
			}else if(fw_visitated == bw_visitated && fw_visitated == false){
				d_colors[v] = 4 * d_pivots[v] + 3;				
			}
				
			if(!get_is_d_eliminated(&d_status[v])){
				*d_stop = false;
			}
		}
		d_write_id_for_pivots[d_colors[v]] = v;
	}
}

__global__ void set_race_winners(int num_nodes, char * d_status, int * d_pivots, int * d_colors, long * d_write_id_for_pivots){
	// Ottenuti i valori della race, si vanno ad impostare i nuovi pivot
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, aggiornata dopo l'esecuzione di update

	
	int v = threadIdx.x + blockIdx.x * blockDim.x;
	if(v < num_nodes) {
		// Se il nodo è stato eliminato, allora il suo pivot è per forza se stesso
		if(get_is_d_eliminated(&d_status[v])){
			d_pivots[v] = v;
		}else{
			d_pivots[v] = d_write_id_for_pivots[d_colors[v]];
			set_is_d_fw_visited(&d_status[d_pivots[v]]);
			set_is_d_bw_visited(&d_status[d_pivots[v]]);
		}
	}
}

__global__ void initialize_pivot(int num_nodes, int * d_pivots, char * d_status) {
	// Scelta iniziale del primo pivot, basandosi sui nodi cancellati inizialmente
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, avente come pivot un nodo non cancellato
	//          fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no. A questo punto l'unico nodo visitato è il solo pivot scelto
	//          bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no. A questo punto l'unico nodo visitato è il solo pivot scelto

	
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		__shared__ int chosen_pivot;
		if(!get_is_d_eliminated(&d_status[v])){
			chosen_pivot = v;
		}

		// Sincronizziamo qui i thread per inizializzare questi array: lanciare un altro thread
		// solo per inizializzare gli array potrebbe risultare più pesante che farlo qui
		__syncthreads();

		d_pivots[v] = chosen_pivot;
		set_is_d_fw_visited(&d_status[d_pivots[v]]);
		set_is_d_bw_visited(&d_status[d_pivots[v]]);
	}
}

void update(int num_nodes, int * d_pivots, char * d_status, long * d_write_id_for_pivots, int * d_colors, bool * stop, const int n_blocks, const int t_per_blocks) {
	// Esegue l'update dei valori del pivot facendo una race
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, aggiornata dopo l'esecuzione di update


	bool *d_stop;

	HANDLE_ERROR(cudaMalloc((void**)&d_stop, sizeof(bool)));
	
	HANDLE_ERROR(cudaMemset(d_write_id_for_pivots, -1, 4 * num_nodes * sizeof(long)));
	HANDLE_ERROR(cudaMemset(d_stop, true, sizeof(bool)));
	
	// Dai paper:
	// These subgraphs are 
	// 		1) the strongly connected component with the pivot;
	// 		2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	// 		3) the subgraph given by vertices in the backward closure but not in the forward closure;
	// 		4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	// The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	// they are recursively processed in parallel with the same algorithm
	
	set_colors<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots, d_stop);
	
	HANDLE_ERROR(cudaMemcpy(stop, d_stop, sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(d_stop));

	// Setto i valori dei pivot che hanno vinto la race
	// Se sono stati eliminati, allora setta il valore dello stesso nodo 
	set_race_winners<<<n_blocks, t_per_blocks>>>(num_nodes, d_status, d_pivots, d_colors, d_write_id_for_pivots);
}

__global__ void trim_u_kernel(int num_nodes, int * d_nodes, int * d_adjacency_list, int * d_pivots, char * d_status, int * d_is_scc){
	// Setta i pivot delle SCC uguale a -1 se questi ricevono archi da nodi u
	// param: 	pivots = 	Lista che per ogni 'v' dice il valore del pivot della SCC
	// 			is_scc =	Lista copia di pivots
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1
	
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		if(get_is_d_u(&d_status[v])){
			for(int u = d_nodes[v]; u < d_nodes[v+1]; ++u) {
				if(d_pivots[v] != d_pivots[d_adjacency_list[u]]) {
					d_is_scc[d_pivots[d_adjacency_list[u]]] = -1;
				}
			}
		}
	}

}

__global__ void trim_u_propagation(int num_nodes, int * d_pivots, int * d_is_scc) {
	// Se alcuni pivot sono settati a -1, per la cancellazione dovuta a collegamenti con nodi u, 
	// propaga la cancellazione agli altri membri della SCC
	// param: 	pivots = 	Lista contenente i pivot delle SCC
	// 			is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot e gli altri nodi delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1

	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes)
		d_is_scc[v] = d_is_scc[d_pivots[v]];
}

__global__ void calculate_more_than_one(int num_nodes, int * d_more_than_one_dev, int * is_scc_dev) {
	// Trova il numero di elementi nella SCC
	// @param: is_scc =	Lista contenente i pivot delle SCC, però i pivot e gli altri nodi delle SCC 
	// 					che ricevono archi da nodi u sono settati a -1
	// @return:	more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 								Se 'v' è pivot: 	more_than_one[v] = numero di elementi nella sua SCC,
	// 								Se 'v' non è pivot:	more_than_one[v] = 1

	int u = threadIdx.x + blockIdx.x * blockDim.x;

	if (u < num_nodes){
		if(is_scc_dev[u] != -1){
			// atomicAdd può essere migliorato -> Simile al problema dell'istogramma
			atomicAdd(&d_more_than_one_dev[is_scc_dev[u]], 1);
		}
	}
}

__global__ void is_scc_adjust(int num_nodes, int * more_than_one_dev, int * is_scc_dev) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// @param: more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 							Se 'v' è pivot: 								more_than_one[v] = numero di elementi nella sua SCC,
	// 							Se 'v' non è pivot, ma fa parte di una SCC:		more_than_one[v] = 0
	// 							Se 'v' non è pivot e non fa parte di una SCC:	more_than_one[v] = 0
	// @return: is_scc =	Lista che per ogni nodo 'v' dice se questo fa parte di una SCC.
	// 						Se fa parte di una SCC: 	is_scc[v] = valore del pivot,
	// 						Se non fa parte di una SCC:	is_scc[v] = -1

	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if(more_than_one_dev[v] == 1)
			is_scc_dev[v] = -1;
	}
}

int count_distinct_scc(int is_scc[], int num_nodes){
	// Restituisce il numero di SCC valide presenti nell'array is_scc
	// Questa funzione non viene parallelizzata poiché utilizzata solamente per verificare la correttezza del risultato
	// @param:  is_scc 	= 	Lista contenente le SCC valide trovate
	// @return: res    	=	Numero di SCC valide diverse

    int res = 0;
 
    // Per tutti gli elementi dell'array
    for (int i = 1; i < num_nodes; i++) {
        int j = 0;
        for (j = 0; j < i; j++)
            if (is_scc[i] == is_scc[j])
                break;
 
        // Se non è già stato contato, contalo
        if (i == j)
            res++;
    }
    return res;
}

int main(int argc, char ** argv) {
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	int num_nodes, num_edges;
    int * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose, * is_scc;
	char * status;

    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);

	if(DEBUG_MAIN){
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_MAIN);
		for (int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_MAIN);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_MAIN);
		for (int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_MAIN);
	}

	const int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const int NUMBER_OF_BLOCKS = num_nodes / THREADS_PER_BLOCK + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);

	// Dichiarazioni di variabili device
	int * d_is_scc, * d_more_than_one, * d_nodes, * d_adjacency_list, * d_nodes_transpose, * d_adjacency_list_transpose, * d_pivots, * d_colors;
	char * d_status;
	long * d_write_id_for_pivots;

	HANDLE_ERROR(cudaMalloc((void**)&d_nodes, (num_nodes+1) * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_nodes_transpose, (num_nodes+1) * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_adjacency_list, num_edges * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_adjacency_list_transpose, num_edges * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_pivots, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_colors, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_status, num_nodes * sizeof(char)));

	HANDLE_ERROR(cudaMalloc((void**)&d_write_id_for_pivots, 4 * num_nodes * sizeof(long)));

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpy(d_nodes, nodes, (num_nodes+1) * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_adjacency_list, adjacency_list, num_edges * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(int), cudaMemcpyHostToDevice));

	HANDLE_ERROR(cudaMemcpy(d_status, status, num_nodes * sizeof(char), cudaMemcpyHostToDevice));

	// Inizializzazione e copia delle funzioni device che verranno passate tramite parametro.
	// Utilizzando le funzioni in questo modo, anche se apparentemente verboso, permette di ottenere meno codice duplicato:
	// infatti, se non fosse per queste variabili, si sarebbe dovuto duplicare l'f_kernel e il reach per averne uno per la forward e uno per la backward.
	get_status h_get_fw_visited, h_get_bw_visited, h_get_fw_expanded, h_get_bw_expanded;
	set_status h_set_fw_visited, h_set_bw_visited, h_set_fw_expanded, h_set_bw_expanded;

	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_fw_visited, dev_get_fw_visited, sizeof(get_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_bw_visited, dev_get_bw_visited, sizeof(get_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_fw_expanded, dev_get_fw_expanded, sizeof(get_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_get_bw_expanded, dev_get_bw_expanded, sizeof(get_status)));

	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_fw_visited, dev_set_fw_visited, sizeof(set_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_bw_visited, dev_set_bw_visited, sizeof(set_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_fw_expanded, dev_set_fw_expanded, sizeof(set_status)));
	HANDLE_ERROR(cudaMemcpyFromSymbol(&h_set_bw_expanded, dev_set_fw_expanded, sizeof(set_status)));
	
	// Primo trimming per eliminare i nodi che, dopo la cancellazione dei nodi non in U,
	// non avevano più out-degree e in-degree diverso da 0
	trimming(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, THREADS_PER_BLOCK, NUMBER_OF_BLOCKS);
	
	// Si fanno competere i thread per scelgliere un nodo che farà da pivot, a patto che quest'ultimo sia non eliminato
	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_status);
	
    bool stop = false;
	
	// Si ripete il ciclo fino a quando tutti i nodi vengono eliminati
    while (!stop){
		// Forward reach
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, h_get_fw_visited, h_get_fw_expanded, h_set_fw_visited, h_set_fw_expanded, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
		
		// Backward reach
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, d_nodes_transpose, d_adjacency_list_transpose, d_pivots, d_status, h_get_bw_visited, h_get_bw_expanded, h_set_bw_visited, h_set_bw_expanded, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);

		// Trimming per eliminare ulteriori nodi che non hanno più out-degree e in-degree diversi da 0
		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, d_nodes, d_nodes_transpose, d_adjacency_list, d_adjacency_list_transpose, d_status, THREADS_PER_BLOCK, NUMBER_OF_BLOCKS);

		// Update dei pivot
		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, d_pivots, d_status, d_write_id_for_pivots, d_colors, &stop, NUMBER_OF_BLOCKS, THREADS_PER_BLOCK);
    }
	
	//Disallocamento della memoria iniziale
	HANDLE_ERROR(cudaFree(d_write_id_for_pivots));
	
	// Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1
	
	// Per iniziare le assegnamo gli stessi valori di pivots, che verranno modificati in seguito
	is_scc = (int*) malloc(num_nodes * sizeof(int));
	
	// Allochiamo more_than_one, che per ogni nodo che fa da pivot viene assegnato un contatore, il quale conta quante volte appare tale pivot
	// Se appare solo una volta, allora il nodo non fa parte di nessuna SCC
	HANDLE_ERROR(cudaMalloc((void**)&d_is_scc, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&d_more_than_one, num_nodes * sizeof(int)));
	
	HANDLE_ERROR(cudaMemcpy(d_is_scc, d_pivots, num_nodes * sizeof(int), cudaMemcpyDeviceToDevice));
	HANDLE_ERROR(cudaMemset(d_more_than_one, 0, num_nodes * sizeof(int)));

	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_nodes, d_adjacency_list, d_pivots, d_status, d_is_scc);
	
	HANDLE_ERROR(cudaFree(d_adjacency_list_transpose));
	HANDLE_ERROR(cudaFree(d_adjacency_list));
	HANDLE_ERROR(cudaFree(d_nodes_transpose));
	HANDLE_ERROR(cudaFree(d_nodes));
	
	trim_u_propagation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_pivots, d_is_scc);

	HANDLE_ERROR(cudaFree(d_pivots));

	calculate_more_than_one<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_more_than_one, d_is_scc);
	is_scc_adjust<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, d_more_than_one, d_is_scc);
	
	HANDLE_ERROR(cudaFree(d_more_than_one));
	HANDLE_ERROR(cudaFree(d_status));

	HANDLE_ERROR(cudaMemcpy(is_scc, d_is_scc, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(d_is_scc));

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("is_scc[" + to_string(i) + "] = ", is_scc[i], false);

	DEBUG_MSG("Number of SCCs found: ", count_distinct_scc(is_scc, num_nodes), DEBUG_FINAL);
}