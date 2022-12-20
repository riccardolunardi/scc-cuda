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

__global__ void f_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status, bool * d_stop, bool (*get_visited)(char *), bool (*get_expanded)(char *), void (*set_visited)(char *), void (*set_expanded)(char *)){
	// Esecuzione di un thread della chiusura in avanti/indietro
	// @param:	pivots			=	Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno
	// 			is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return 	is_visited		=	Lista che per ogni 'v' dice se è stato visitato dalla reach o meno, aggiornata dopo l'esecuzione del trimming
	// 			is_expanded		=	Lista che per ogni 'v' dice se sono stato visitati i figli diretti o meno, aggiornata dopo l'esecuzione del trimming
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

    // Per ogni nodo
	if(v < num_nodes) {
        // Si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		char node_status = d_status[v];
		if(!get_is_d_eliminated(&node_status) && get_visited(&node_status) && !get_expanded(&node_status)) {
            // Si segna come espanso
			set_expanded(&d_status[v]);

            // Per ogni nodo a cui punta
			// Nonostrante il for crei molta divergenza, nei grandi grafi è preferibile utilizzare
			// questo metodo che, ad esempio, il dynamic programming.
			for(unsigned int u = d_nodes[v]; u < d_nodes[v + 1]; u++) {	
				unsigned int dst = d_adjacency_list[u];
				char dst_status = d_status[dst];

                // Si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!get_is_d_eliminated(&dst_status) && !get_visited(&dst_status)) {
                    // Setta il nodo puntato a visitato
					set_visited(&d_status[dst]);
                    // Permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					*d_stop = false;
				}
			}
		}
	}
}

__global__ void bitwise_or_kernel(const unsigned int num_nodes, char * d_status_res, char * d_bw_status){
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		d_status_res[v] |= d_bw_status[v];
	}
}

__global__ void trimming_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_nodes_transpose, unsigned int * d_adjacency_list,  unsigned int * d_adjacency_list_transpose, char * d_status, bool * d_stop){
	// Esegue un'eliminazione di nodi con out-degree o in-degree uguale a 0, senza contare i nodi eliminati
	// @param:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// @return:	is_eliminated	=	Lista che per ogni 'v' dice se il nodo è stato eliminato o no, aggiornata dopo l'esecuzione del trimming
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		if(!get_is_d_eliminated(&d_status[v])){
			// Se questo valore non verrà cambiato, allora il nodo verrà cancellato
			bool elim = true;
			bool forward = false;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0, tra i soli nodi non eliminati, allora non va eliminato
			unsigned int src = d_nodes[v];
			unsigned int dst = d_nodes[v + 1];
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

			if(elim){
				set_is_d_eliminated(&d_status[v]);
				*d_stop = false;
			}
		}
	}
}

__global__ void set_colors(unsigned int const num_nodes, char * d_status, unsigned int * d_pivots, unsigned int * d_colors, unsigned long * d_write_id_for_pivots, bool * d_stop){
	// Esegue l'update dei valori del pivot facendo una race, scrivendo il "colore" di una serie di pivot in array simultaneamente
	// @param:	pivots						= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated				= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: d_write_id_for_pivots		= Lista che conterrà, nelle posizione identificate dai colori appena calcolati, i nuovi pivot da assegnare
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		unsigned int new_color;
		char src = d_status[v];

		if(!get_is_d_eliminated(&src)){
			unsigned int pivot_node = d_pivots[v];
			const bool fw_visitated = get_is_d_fw_visited(&src);
			const bool bw_visitated = get_is_d_bw_visited(&src);
			
			if(fw_visitated == bw_visitated && fw_visitated == true){
				new_color = pivot_node << 2;
			} else {
				//printf("Nodo %d ha fw_visitated = %d e bw_visitated = %d\n", v, fw_visitated, bw_visitated);
				*d_stop = false;

				if(fw_visitated != bw_visitated && fw_visitated == true){
					new_color = (pivot_node << 2) + 1;
				}else if(fw_visitated != bw_visitated && fw_visitated == false){
					new_color = (pivot_node << 2) + 2;
				}else if(fw_visitated == bw_visitated && fw_visitated == false){
					new_color = (pivot_node << 2) + 3;				
				}
			}

			d_write_id_for_pivots[new_color] = v;
			d_colors[v] = new_color;
		}

		/* __syncthreads();

		// Se il nodo è stato eliminato, allora il suo pivot è per forza se stesso
		if(get_is_d_eliminated(&src)){
			if(!get_is_d_scc(&src)){
				d_pivots[v] = v;
			}
		}else{
			d_pivots[v] = d_write_id_for_pivots[new_color];
			set_is_d_bw_fw_visited(&d_status[d_pivots[v]]);		
			if(new_color % 4 == 0){
				set_is_d_eliminated(&d_status[v]);
				set_is_d_scc(&d_status[v]);
			}
		} */
		

	}
}

__global__ void set_new_pivots(unsigned int const num_nodes, char * d_status, unsigned int * d_pivots, unsigned int * d_colors, unsigned long * d_write_id_for_pivots){
	// Esegue l'update dei valori del pivot facendo una race, scrivendo il "colore" di una serie di pivot in array simultaneamente
	// @param:	pivots						= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated				= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited				= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: d_write_id_for_pivots		= Lista che conterrà, nelle posizione identificate dai colori appena calcolati, i nuovi pivot da assegnare
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		char src = d_status[v];
		const unsigned int new_color = d_colors[v];

		// Se il nodo è stato eliminato, allora il suo pivot è per forza se stesso
		if(get_is_d_eliminated(&src)){
			if(!get_is_d_scc(&src)){
				d_pivots[v] = v;
			}
		}else{
			d_pivots[v] = d_write_id_for_pivots[new_color];
			set_is_d_bw_fw_visited(&d_status[d_pivots[v]]);		
			if(new_color % 4 == 0){
				set_is_d_scc_is_eliminated(&d_status[d_pivots[v]]);
				set_is_d_scc_is_eliminated(&d_status[v]);
			}
		}
		

	}
}

__global__ void initialize_pivot(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status) {
	// Scelta iniziale del primo pivot, basandosi sui nodi cancellati inizialmente
	// @param:	pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene
	// 			is_eliminated	= Lista che per ogni 'v' dice se il nodo è stato eliminato o no
	// 			fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no
	// 			bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no
	// @return: pivots			= Lista che contiene, per ogni 'v', il valore del pivot della SCC a cui tale nodo 'v' appartiene, avente come pivot un nodo non cancellato
	//          fw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la forward reach partendo dai pivots o no. A questo punto l'unico nodo visitato è il solo pivot scelto
	//          bw_is_visited	= Lista che per ogni 'v' dice se il nodo è stato visitato con la backward reach partendo dai pivots o no. A questo punto l'unico nodo visitato è il solo pivot scelto

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		//__shared__ unsigned int chosen_pivot;
		//chosen_pivot = 2000000000;

		if(!get_is_d_eliminated(&d_status[v])){
			d_pivots[0] = v;
		}
	}
}

__global__ void set_initialize_pivot(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status) {
	// Sincronizziamo qui i thread del blocco per inizializzare questi array: lanciare un altro thread
	// solo per inizializzare gli array potrebbe risultare più pesante che farlo qui
	/* __syncthreads();

	if(threadIdx.x == 0){
		d_pivots[0] = chosen_pivot;
	} */
	
	//__syncthreads();*/
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		d_pivots[v] = d_pivots[0];
		set_is_d_bw_fw_visited(&d_status[d_pivots[0]]);
	}
}

__global__ void trim_u_kernel(unsigned int const num_nodes, unsigned int * d_nodes, unsigned int * d_adjacency_list, unsigned int * d_pivots, char * d_status){
	// Setta i pivot delle SCC uguale a -1 se questi ricevono archi da nodi u
	// param: 	pivots = 	Lista che per ogni 'v' dice il valore del pivot della SCC
	// 			is_scc =	Lista copia di pivots
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1
	
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		unsigned int v_pivot = d_pivots[v];

		if(get_is_d_u(&d_status[v])){
			unsigned int start = d_nodes[v];
			unsigned int lst = d_nodes[v+1];

			for(unsigned int u = start; u < lst; ++u) {
				unsigned int dst = d_adjacency_list[u];
				unsigned int dst_pivot = d_pivots[dst];

				if(v_pivot != dst_pivot) {
					set_not_is_d_scc(&d_status[dst_pivot]);
				}
			}
		}
	}

}

__global__ void trim_u_propagation(unsigned int const num_nodes, unsigned int * d_pivots, char * d_status, bool * is_scc) {
	// Se alcuni pivot sono settati a -1, per la cancellazione dovuta a collegamenti con nodi u, 
	// propaga la cancellazione agli altri membri della SCC
	// param: 	pivots = 	Lista contenente i pivot delle SCC
	// 			is_scc =	Lista contenente i pivot delle SCC, però i pivot delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1
	// @return:	is_scc =	Lista contenente i pivot delle SCC, però i pivot e gli altri nodi delle SCC 
	// 						che ricevono archi da nodi u sono settati a -1

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

__global__ void is_scc_adjust(unsigned int const num_nodes, int unsigned * d_pivots, char * d_status) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// @param: more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 							Se 'v' è pivot: 								more_than_one[v] = numero di elementi nella sua SCC,
	// 							Se 'v' non è pivot, ma fa parte di una SCC:		more_than_one[v] = 0
	// 							Se 'v' non è pivot e non fa parte di una SCC:	more_than_one[v] = 0
	// @return: is_scc =	Lista che per ogni nodo 'v' dice se questo fa parte di una SCC.
	// 						Se fa parte di una SCC: 	is_scc[v] = valore del pivot,
	// 						Se non fa parte di una SCC:	is_scc[v] = -1

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){

		if(d_pivots[v] == v)
			set_not_is_d_scc(&d_status[v]);

		/* __syncthreads();

		if (!get_is_d_scc(&d_status[v]))
			set_not_is_d_scc(&d_status[d_pivots[v]]); */
	}
}

__global__ void is_scc_adjust_prop(unsigned int const num_nodes, int unsigned * d_pivots, char * d_status) {
	// Restituisce una lista che dice se il nodo 'v' fa parte di una SCC
	// @param: more_than_one = 	Lista che per ogni nodo 'v' dice se questo è un pivot.
	// 							Se 'v' è pivot: 								more_than_one[v] = numero di elementi nella sua SCC,
	// 							Se 'v' non è pivot, ma fa parte di una SCC:		more_than_one[v] = 0
	// 							Se 'v' non è pivot e non fa parte di una SCC:	more_than_one[v] = 0
	// @return: is_scc =	Lista che per ogni nodo 'v' dice se questo fa parte di una SCC.
	// 						Se fa parte di una SCC: 	is_scc[v] = valore del pivot,
	// 						Se non fa parte di una SCC:	is_scc[v] = -1

	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		if (!get_is_d_scc(&d_status[v]))
			set_not_is_d_scc(&d_status[d_pivots[v]]);
	}
}

__global__ void eliminate_trivial_scc(unsigned int const t_p_b, unsigned int const num_nodes, int unsigned * d_pivots, bool * d_is_scc) {
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	extern __shared__ bool s_scc_pivots[];
	bool *s_is_scc = s_scc_pivots;
	unsigned int *s_pivots = (unsigned int*)&s_scc_pivots[t_p_b];

	if (v < num_nodes){
		s_is_scc[threadIdx.x] = d_is_scc[v];
		s_pivots[threadIdx.x] = d_pivots[v];

		__syncthreads();

		if(s_pivots[threadIdx.x] == v)
			s_is_scc[threadIdx.x] = false;

		__syncthreads();

		d_is_scc[v] = s_is_scc[threadIdx.x];
		d_pivots[v] = s_pivots[threadIdx.x];
	}
}

__global__ void convert_int_array_to_bool(unsigned int const num_nodes, int * d_is_scc, bool * d_is_scc_final) {
	unsigned int const v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		d_is_scc_final[v] = d_is_scc[v] != -1;
	}
}

__global__ void is_scc_adjust_v1(int num_nodes, int * more_than_one_dev, int * is_scc_dev) {
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

unsigned count_distinct_scc(unsigned n, unsigned int pivots[], char status[]){
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

int count_distinct_scc_v1(int is_scc[], int num_nodes){
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

__global__ void at_least_one_scc(const unsigned int num_nodes, bool * d_is_scc){
	unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < num_nodes){
		if (d_is_scc[i] > 0){
			d_is_scc[0] = true;
		}
	}
}

bool is_there_an_scc(const unsigned int NUMBER_OF_BLOCKS, const unsigned int thread_per_block, const unsigned int num_nodes, bool * d_is_scc){
	bool final_result = false;

	at_least_one_scc<<<NUMBER_OF_BLOCKS, thread_per_block>>>(num_nodes, d_is_scc);
	HANDLE_ERROR(cudaMemcpy(&final_result, (bool*)d_is_scc, sizeof(bool), cudaMemcpyDeviceToHost));

	return final_result;
}

#endif