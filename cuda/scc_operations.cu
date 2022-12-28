#ifndef SCC_OPERATIONS
#define SCC_OPERATIONS

#include <cuda.h>
#include <set>
#include "../utils/is_checked.cu"
#include "../utils/is_checked.cpp"

static void handle_error(cudaError_t err, const char *file, int line ) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
		exit( EXIT_FAILURE );
	}
}
#define HANDLE_ERROR( err ) (handle_error( err, __FILE__, __LINE__ ))

inline __host__ __device__ void operator|=(char4 &a, char4 b){
	// Definizione dell'operatore di assegnamente bit a bit OR per le variabili di tipo char4
	// Il tipo char4 è un tipo di dato che contiene 4 variabili di tipo char, definito in cuda_runtime.h
	// La funzione è definita inline per permettere al compilatore di sostituire il codice della funzione
    a.x |= b.x;
    a.y |= b.y;
    a.z |= b.z;
    a.w |= b.w;
}

inline __host__ __device__ bool operator>(const char4 &a, const int b){
	// Definizione dell'operatore > per le variabili di tipo char4. Se la somma dei valori è maggiore di b, allora ritorna true
	// Il tipo char4 è un tipo di dato che contiene 4 variabili di tipo char, definito in cuda_runtime.h
	// La funzione è definita inline per permettere al compilatore di sostituire il codice della funzione
    return a.x + a.y + a.z + a.w > b;
}

__global__ void f_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status, bool * d_stop, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *)){
	// Esecuzione di un thread della chiusura in avanti/indietro
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

    // Per ogni nodo v
	if(v < num_nodes) {
        // Si controlla se v non è stato eliminato, se è stato visitato e se non è stato espanso
		char node_status = d_status[v];
		if(!get_is_d_eliminated(&node_status) && get_visited(&node_status) && !get_expanded(&node_status)) {
            // Si segna come espanso, per non ricontrollarlo più nelle prossime iterazioni
			set_expanded(&d_status[v]);

            // Per ogni nodo u a cui punta...
			// Nonostrante il for crei molta divergenza, nei grandi grafi è preferibile utilizzare
			// questo metodo che, ad esempio, il dynamic programming.
			for(unsigned int u = d_nodes[v]; u < d_nodes[v + 1]; u++) {	
				unsigned int dst = d_adjacency_list[u];
				char dst_status = d_status[dst];

                // Si controlla se u non è stato eliminato e se non è stato visitato
				if(!get_is_d_eliminated(&dst_status) && !get_visited(&dst_status)) {
                    // Setta il nodo u come visitato
					set_visited(&d_status[dst]);
                    // Si è trovato un altro nodo visitato ancora da espandere, quindi continuo il ciclo reach
					*d_stop = false;
				}
			}
		}
	}
}

__global__ void bitwise_or_kernel(const unsigned int num_nodes, char * d_status_res, char * d_bw_status){
	// Esegue un'operazione di OR bit a bit tra i due array contenenti rispettivamente i risultati della chiusura in avanti e indietro
	// Vista la semplicità dell'operazione, è stata implementata in modo da utilizzare l'accesso vettorizzato
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		const int node_limit = num_nodes >> 2;
		for(int i = v; i < node_limit; i += blockDim.x * gridDim.x) {
			reinterpret_cast<char4*>(d_status_res)[v] |= reinterpret_cast<char4*>(d_bw_status)[v];
		}

		// In un solo thread si controlla se il numero di nodi è un multiplo di 4,
		// In caso contrario si esegue l'operazione di OR singolarmente
		int remainder = num_nodes & 3;
		if (v == node_limit && remainder != 0) {
			while(remainder) {
				int idx = num_nodes - remainder--;
				d_status_res[v] |= d_bw_status[v];
			}
	  }
	}
}

__global__ void trimming_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list,  unsigned int * d_adjacency_list_transpose, char * d_status, bool * d_stop){
	// Esegue un'eliminazione di nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		// Per ogni nodo v
		if(!get_is_d_eliminated(&d_status[v])){
			// Se non è stato eliminato
			bool elim = true;
			bool forward = false;
			
			unsigned int src = d_nodes[v];
			unsigned int dst = d_nodes[v + 1];
			// Se v, contando solo i nodi non eliminati, ha sia in_degree > 0 che out_degree > 0 allora non va eliminato
			for(unsigned int u = src; u < dst; u++){
				if(!get_is_d_eliminated(&d_status[d_adjacency_list[u]])) {
					forward = true;
				}
			}
			if(forward) {
				src = d_nodes_transpose[v];
				dst = d_nodes_transpose[v + 1];

				for(unsigned int u = src; u < dst; u++){
					if(!get_is_d_eliminated(&d_status[d_adjacency_list_transpose[u]])) {
						elim = false;
					}
				}
			}

			// Se elim non è stato modificato, allora il nodo v va eliminato
			if(elim){
				set_is_d_eliminated(&d_status[v]);
				*d_stop = false;
			}
		}
	}
}

__global__ void set_colors(unsigned int const num_nodes, char * d_status, unsigned int * d_pivots, unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * d_stop){
	// Esegue l'update dei valori del pivot facendo una race, scrivendo il "colore" di una serie di pivot in array simultaneamente
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	// Per ogni nodo v
	if(v < num_nodes) {
		unsigned int new_color;
		char src = d_status[v];

		// Se non è stato eliminato
		if(!get_is_d_eliminated(&src)){
			unsigned int pivot_node = d_pivots[v];
			const bool fw_visitated = get_is_d_fw_visited(&src);
			const bool bw_visitated = get_is_d_bw_visited(&src);
			
			// Se fa parte di una SCC, quindi è stato visitato sia in avanti che all'indietro
			if(fw_visitated == bw_visitated && fw_visitated == true){
				new_color = pivot_node << 2;
			} else {
				*d_stop = false;

				// Se è stato visitato solo in avanti
				if(fw_visitated != bw_visitated && fw_visitated == true){
					new_color = (pivot_node << 2) + 1;
				// Se è stato visitato solo all'indietro
				}else if(fw_visitated != bw_visitated && fw_visitated == false){
					new_color = (pivot_node << 2) + 2;
				// Se non è stato visitato né in avanti né all'indietro
				}else if(fw_visitated == bw_visitated && fw_visitated == false){
					new_color = (pivot_node << 2) + 3;				
				}
			}

			// Su questa riga viene effettuata la race
			// Ogni sotto-grafo avrà come pivot l'ultimo nodo che esegue questa riga
			d_write_id_for_pivots[new_color] = v;
			d_colors[v] = new_color;
		}
	}
}

__global__ void set_new_pivots(unsigned int const num_nodes, char * d_status, unsigned int * d_pivots, unsigned int * d_colors, unsigned long * d_write_id_for_pivots){
	// Esegue l'update dei valori del pivot facendo una race, scrivendo il "colore" di una serie di pivot in array simultaneamente
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		if(!get_is_d_eliminated(&d_status[v])){
			// Se non sono stati eliminati, allora setta il valore del pivot uguale al nodo che ha vinto la race
			d_pivots[v] = d_write_id_for_pivots[d_colors[v]];
			set_is_d_bw_fw_visited(&d_status[d_pivots[v]]);		
			// I nodi che fanno parte di una SCC, vengono settati come eliminati e come SCC
			if((d_colors[v] & 3) == 0){
				set_is_d_scc(&d_status[d_pivots[v]]);
				set_is_d_scc(&d_status[v]);
			}
		}
	}
}

__global__ void set_new_eliminated(unsigned int const num_nodes, char * d_status, unsigned int * d_pivots, unsigned int * d_colors, unsigned long * d_write_id_for_pivots){

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		if(get_is_d_eliminated(&d_status[v])){
			if(!get_is_d_scc(&d_status[v])){
				d_pivots[v] = v;
			}
		}
		if(get_is_d_scc(&d_status[v])){
			set_is_d_eliminated(&d_status[v]);
		}
	}

}

__global__ void initialize_pivot(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status) {
	// Scelta iniziale del primo pivot, basandosi sui nodi cancellati inizialmente

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		// Se un nodo non è stato eliminato, allora il suo pivot proverà ad essere scritto nella
		// prima posizione dell'array dei pivot. L'ultimo pivot scritto in questa posizione sarà
		// Il primo pivot iniziale di tutti
		if(!get_is_d_eliminated(&d_status[v])){
			d_pivots[0] = v;
		}
	}
}

__global__ void set_initialize_pivot(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status) {
	// Sincronizziamo qui i thread del blocco per inizializzare questi array: lanciare un altro thread
	// solo per inizializzare gli array potrebbe risultare più pesante che farlo qui
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		const int node_limit = num_nodes >> 2;
		for(int i = v; i < node_limit; i += blockDim.x * gridDim.x) {
			reinterpret_cast<uint4*>(d_pivots)[v] = make_uint4(d_pivots[0], d_pivots[0], d_pivots[0], d_pivots[0]);
		}

		// in only one thread, process final elements (if there are any)
		int remainder = num_nodes & 3;
		if (v==node_limit && remainder!=0) {
			// L'operazione di settaggio di fw e bw visited è stata spostata qui per evitare di farla in ogni thread
			// È infatti sufficiente farla una sola volta, perché a questo punto il pivot è uno solo
			set_is_d_bw_fw_visited(&d_status[d_pivots[v]]);
			
			while(remainder) {
				int idx = num_nodes - remainder--;
				d_pivots[v] = d_pivots[0];
			}
	  	}
	}
}

__global__ void trim_u_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status){
	// Controlla per tutte le SCC se queste ricevono archi da nodi u
	// Se le trova setta i pivot come non facenti parte di una SCC
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		unsigned int v_pivot = d_pivots[v];

		if(get_is_d_u(&d_status[v])){
			unsigned int start = d_nodes[v];
			unsigned int lst = d_nodes[v+1];

			for(unsigned int u = start; u < lst; ++u) {
				unsigned int dst = d_adjacency_list[u];
				unsigned int dst_pivot = d_pivots[dst];

				// Dato un arco (u,v), se u non fa parte della stessa SCC di v e u fa parte di U
				// Allora setto il pivot della SCC di v come non facente parte di una SCC
				if(v_pivot != dst_pivot) {
					set_not_is_d_scc(&d_status[dst_pivot]);
				}
			}
		}
	}

}

__global__ void trim_u_propagation(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status, bool * is_scc) {
	// Se sono presenti pivot non facenti più parte di una SCC, per la cancellazione dovuta a trim_u_kernel, 
	// propaga la cancellazione agli altri nodi della stessa SCC

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if(get_is_d_scc(&d_status[d_pivots[v]])){
			set_is_d_scc(&d_status[v]);
			is_scc[v] = true;
		}else{
			set_not_is_d_scc(&d_status[v]);
			is_scc[v] = false;
		}
	}
}

__global__ void trim_u_propagation_v1(int num_nodes, int * d_pivots, int * d_is_scc) {
	// Se sono presenti pivot non facenti più parte di una SCC, per la cancellazione dovuta a trim_u_kernel, 
	// propaga la cancellazione agli altri nodi della stessa SCC
	// Questa versione è utilizzata per la versione 1 dell'implementazione dell'algoritmo in CUDA

	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes)
		d_is_scc[v] = d_is_scc[d_pivots[v]];
}

__global__ void is_scc_adjust_v1(int num_nodes, int * more_than_one_dev, int * is_scc_dev) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// Questa versione è utilizzata per la versione 1 dell'implementazione dell'algoritmo in CUDA

	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if(more_than_one_dev[v] == 1)
			is_scc_dev[v] = -1;
	}
}

__global__ void is_scc_adjust(unsigned int const num_nodes, int unsigned * d_pivots, char * d_status) {
	// Modifica lo status di un nodo v nel caso quest'ultimo non faccia parte di una SCC

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if(d_pivots[v] == v)
			set_not_is_d_scc(&d_status[v]);
	}
}

__global__ void is_scc_adjust_prop(unsigned int const num_nodes, int unsigned * d_pivots, char * d_status) {
	// Propagazione del risultato di is_scc_adjust

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if (!get_is_d_scc(&d_status[v]))
			set_not_is_d_scc(&d_status[d_pivots[v]]);
	}
}

__global__ void eliminate_trivial_scc(unsigned int const t_p_b, unsigned int const num_nodes, int unsigned * d_pivots, bool * d_is_scc) {
	// Setta tutti i pivot come non facenti parte di una SCC
	// In questa funzione viene utilizzata la memoria shared per risparmiare tempo di accesso alla memoria globale
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	// Definizione della shared memory
	// Visto che nella chiamata del kernel è obbligatorio dichiarare tutta insieme la memoria shared, bisogna creare un unico array.
	// L'unico array verrà diviso in due parti, una per i pivot e una per i booleani che indicano se un nodo fa parte di una SCC
	// Per poter rendere il codice più leggibile, si definiscono due puntatori che puntano all'inizio di ciascuna parte dell'array
	extern __shared__ bool s_scc_pivots[];
	bool *s_is_scc = s_scc_pivots;
	unsigned int *s_pivots = (unsigned int*)&s_scc_pivots[t_p_b];

	if (v < num_nodes){
		// Caricamento dei dati nella shared memory
		s_is_scc[threadIdx.x] = d_is_scc[v];
		s_pivots[threadIdx.x] = d_pivots[v];

		// Barriera di sincronizzazione (per aspettare che tutti i thread abbiano caricato i dati in memoria)
		__syncthreads();

		// Se il pivot di v è uguale a v, allora v non fa parte di una SCC (codice effettivo del kernel)
		if(s_pivots[threadIdx.x] == v)
			s_is_scc[threadIdx.x] = false;

		// Barriera di sincronizzazione (per poter copiare i dati dalla shared memory alla memoria globale)
		__syncthreads();

		d_is_scc[v] = s_is_scc[threadIdx.x];
		d_pivots[v] = s_pivots[threadIdx.x];
	}
}

__global__ void convert_int_array_to_bool(unsigned int const num_nodes, int * d_is_scc, bool * d_is_scc_final) {
	// Funzione che converte un array di interi in un array di booleani
	// La funzione non è ottimizzata perché non strettamente necessario, ma è stata utilizzata per fare debug
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		d_is_scc_final[v] = d_is_scc[v] != -1;
	}
}

/* int count_distinct_scc_v1(int is_scc[], int num_nodes){
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
} */

set<int> count_distinct_scc_v1(int is_scc[], int n){
	// Conta quanti elementi distinti ci sono in un array

	set<int> s;

	// Aggiungo un elemento al set se fa parte della SCC
	// set non permette elementi ripetuti, quindi ogni pivot comparirà una volta sola

	for(int i=0; i<n; i++) {
		if(is_scc[i] != -1) {
        	s.insert(is_scc[i]);
		}
    }
	
    return s;
}

set<unsigned> count_distinct_scc(unsigned n, unsigned int pivots[], char status[]){
	// Conta quanti elementi distinti ci sono in un array

	set<unsigned> s;

	// Aggiungo un elemento al set se fa parte della SCC
	// set non permette elementi ripetuti, quindi ogni pivot comparirà una volta sola

	for(int i=0; i<n; i++) {
		if(get_is_scc(status[i])) {
        	s.insert(pivots[i]);
		}
    }
	
    return s;
}

__global__ void calculate_more_than_one(int num_nodes, int * d_more_than_one_dev, int * is_scc_dev) {
	// Trova il numero di elementi nella SCC
	// Funzione non ottimizzata, ma utilizzata per fare debug nella versione naive di CUDA

	int u = threadIdx.x + blockIdx.x * blockDim.x;

	if (u < num_nodes){
		if(is_scc_dev[u] != -1){
			atomicAdd(&d_more_than_one_dev[is_scc_dev[u]], 1);
		}
	}
}

__global__ void at_least_one_scc(const unsigned int num_nodes, bool * d_is_scc){
	// Viene eseguita una race nella prima cella di memoria per verificare se esiste almeno una SCC
	// Vista la semplicità dell'operazione, è stata implementata in modo da utilizzare l'accesso vettorizzato
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		const int node_limit = num_nodes >> 2;
		for(int i = v; i < node_limit; i += blockDim.x * gridDim.x) {
			if(reinterpret_cast<char4*>(d_is_scc)[i] > 0){
				d_is_scc[0] = true;
			}
		}

		// In un solo thread si controlla se il numero di nodi è un multiplo di 4,
		// In caso contrario si esegue l'operazione di OR singolarmente
		int remainder = num_nodes & 3;
		if (v == node_limit && remainder != 0) {
			while(remainder) {
				int idx = num_nodes - remainder--;
				if (d_is_scc[idx] > 0) {
					d_is_scc[0] = true;
				}
			}
	  }
	}
}

bool is_there_an_scc(const unsigned int NUMBER_OF_BLOCKS, const unsigned int thread_per_block, const unsigned int num_nodes, bool * d_is_scc){
	// Funzione che controlla se esiste almeno una SCC e salva i risultato in una variabile booleana
	bool final_result = false;

	at_least_one_scc<<<NUMBER_OF_BLOCKS, thread_per_block>>>(num_nodes, d_is_scc);
	HANDLE_ERROR(cudaMemcpy(&final_result, (bool*)d_is_scc, sizeof(bool), cudaMemcpyDeviceToHost));

	return final_result;
}

#endif