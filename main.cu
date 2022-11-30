#include "utils.cpp"
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

__global__ void f_kernel(int num_nodes, int num_edges, int * dev_nodes, int * dev_adjacency_list, int * dev_pivots, bool * dev_is_visited, bool * dev_is_eliminated, bool * dev_is_expanded, bool * dev_stop){
    /* for (int i = 0; i < num_nodes; i++){
        DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_F_KERNEL);
    } */
	int v = threadIdx.x + blockIdx.x * blockDim.x;

    // per ogni nodo
	if(v < num_nodes) {
		/* printf("Checking %d\n", v);
        printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
		printf("is_visited[%d] -> %d\n", v, dev_is_visited[v]);
		printf("is_expanded[%d] -> %d\n", v, dev_is_expanded[v]); */
		
        // si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!dev_is_eliminated[v] && dev_is_visited[v] && !dev_is_expanded[v]) {
            // si segna come espanso
			dev_is_expanded[v] = true;
			// printf(" u va da %d a %d\n", dev_nodes[v], dev_nodes[v+1]);
            // per ogni nodo a cui punta
			for(int u = dev_nodes[v]; u < dev_nodes[v + 1]; u++) {	
				/* if(v==138){
					printf("  u -> %d\n", u);
					printf("  dev_adjacency_list[%d] -> %d\n", u, dev_adjacency_list[u]);
				} */
				/* printf("  Nodo %d connesso a nodo %d\n", v, dev_adjacency_list[u]);	
				*/
				int dst = dev_adjacency_list[u];
				/* if(dst==128){
					printf("  dev_is_eliminated[%d] -> %d\n", dst, dev_is_eliminated[dst]);
					printf("  dev_is_visited[%d] -> %d\n", dst, dev_is_visited[dst]);
					printf("  dev_pivots[%d] == dev_pivots[%d] -> %d == %d\n", v, dst, dev_pivots[v], dev_pivots[dst]);
				} */

                // si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				if(!dev_is_eliminated[dst] && !dev_is_visited[dst] && dev_pivots[v] == dev_pivots[dst]) {
                    // setta il nodo puntato a visitato
					//printf("   dev_is_visited[%d] -> %d\n", dst, dev_is_visited[dst]);
					dev_is_visited[dst] = true;
					//printf("   dev_is_visited[%d] -> %d\n", dst, dev_is_visited[dst]);
                    // permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					*dev_stop = false;
					//printf("   %d\n", *dev_stop);
				}
			}
		}
	}
}

__global__ void set_pivots_visited(int num_nodes, bool * dev_is_visited, int * dev_pivots){
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes){
		dev_is_visited[dev_pivots[v]] = true;
	}
}

void reach(int num_nodes, int num_edges, int * dev_nodes, int * dev_adjacency_list, int * dev_pivots, bool * is_visited, bool * dev_is_eliminated, bool * is_expanded) {
    int NOB = num_nodes / 1024 + (num_nodes % 1024 == 0 ? 0 : 1);

	bool *dev_is_visited, *dev_is_expanded;
	bool stop, *dev_stop;
	stop = false;

	HANDLE_ERROR(cudaMalloc((void**)&dev_is_visited, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_expanded, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));

	HANDLE_ERROR(cudaMemcpy(dev_is_visited, is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_expanded, is_expanded, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));

	// Tutti i pivot vengono segnati come visitati
	set_pivots_visited<<<NOB, 1024>>>(num_nodes, dev_is_visited, dev_pivots);
	
    // si effettua la chiusura in avanti/indietro
    while(!stop) {
		stop = true;
		HANDLE_ERROR(cudaMemcpy(dev_stop, &stop, sizeof(bool), cudaMemcpyHostToDevice));
		
        f_kernel<<<NOB, 1024>>>(num_nodes, num_edges, dev_nodes, dev_adjacency_list, dev_pivots, dev_is_visited, dev_is_eliminated, dev_is_expanded, dev_stop);	

		HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }
	
	HANDLE_ERROR(cudaMemcpy(is_visited, dev_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(is_expanded, dev_is_expanded, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));

	for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_visited[" << i << "] -> ", is_visited[i], DEBUG_REACH);
	}

	HANDLE_ERROR(cudaFree(dev_is_visited));
	HANDLE_ERROR(cudaFree(dev_is_expanded));
	HANDLE_ERROR(cudaFree(dev_stop));

}

__global__ void trimming_kernel(int num_nodes, int * dev_nodes, int * dev_nodes_transpose, int * dev_adjacency_list,  int * dev_adjacency_list_transpose, bool * dev_is_eliminated, bool * dev_stop){
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		//printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
		if(!dev_is_eliminated[v]){
			bool elim = true;

			bool forward = false;
			bool backward = false;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0, tra i soli nodi non eliminati, allora non va eliminato
			for(int u = dev_nodes[v]; u < dev_nodes[v+1]; u++){
				if(!dev_is_eliminated[dev_adjacency_list[u]]) {
					forward = true;
				}
			}
			if(forward) {
				for(int u = dev_nodes_transpose[v]; u < dev_nodes_transpose[v+1]; u++){
					if(!dev_is_eliminated[dev_adjacency_list_transpose[u]]) {
						backward = true;
					}
				}
			}
			if(backward) {
				elim = false;
			}

			if(elim){
				dev_is_eliminated[v] = true;
				//printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
				*dev_stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int num_edges, int * dev_nodes, int * dev_nodes_transpose, int * dev_adjacency_list, int * dev_adjacency_list_transpose, bool * dev_is_eliminated, int NUMBER_OF_BLOCKS, int THREADS_PER_BLOCK) {
	bool stop, *dev_stop;
	stop = false;


	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));
	/* for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	} */

    while(!stop) {
		HANDLE_ERROR(cudaMemset(dev_stop, true, sizeof(bool)));
        trimming_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_nodes, dev_nodes_transpose, dev_adjacency_list, dev_adjacency_list_transpose, dev_is_eliminated, dev_stop);
		HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
    }

	HANDLE_ERROR(cudaFree(dev_stop));

	/* for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	} */
}

__global__ void set_colors(int num_nodes, bool * dev_fw_is_visited, bool * dev_bw_is_visited, int * dev_pivots, int * dev_colors, bool * dev_is_eliminated, long * dev_write_id_for_pivots, bool * dev_stop){
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes) {
		// Questo primo caso non ha senso di esistere, perché possiamo lasciargli il valore precedente, tanto cambiaeranno tutti gli altri
		// in realtà ha senso per conservare il valore del pivot, se poi si scopre che una volta diventato SCC il suo valore nel vettore pivots, allora il primo caso si può cancellare e moltiplicare per 3
		if(dev_is_eliminated[v]){
			dev_pivots[v] = v;
		} 
		
		if(dev_fw_is_visited[v] == dev_bw_is_visited[v] && dev_fw_is_visited[v] == true){
			dev_colors[v] = 4 * dev_pivots[v];
		} else {
			if(dev_fw_is_visited[v] != dev_bw_is_visited[v] && dev_fw_is_visited[v] == true){
				dev_colors[v] = 4 * dev_pivots[v] + 1;
			}else if(dev_fw_is_visited[v] != dev_bw_is_visited[v] && dev_fw_is_visited[v] == false){
				dev_colors[v] = 4 * dev_pivots[v] + 2;
			}else if(dev_fw_is_visited[v] == dev_bw_is_visited[v] && dev_fw_is_visited[v] == false){
				dev_colors[v] = 4 * dev_pivots[v] + 3;				
			}
				
			if(!dev_is_eliminated[v]){
				*dev_stop = false;
				//printf("%d -> non eliminato, ma non visitato da fw e bw\n", v);
			}
		}
		dev_write_id_for_pivots[dev_colors[v]] = v;
	}
}

__global__ void set_race_winners(int num_nodes, bool * dev_is_eliminated, int * dev_pivots, int * dev_colors, long * dev_write_id_for_pivots){
	int v = threadIdx.x + blockIdx.x * blockDim.x;
	if(v < num_nodes) {
		if(dev_is_eliminated[v]){
			dev_pivots[v] = v;
		}else{
			dev_pivots[v] = dev_write_id_for_pivots[dev_colors[v]];
		}
	}
}

__global__ void initialize_pivot(int num_nodes, bool * dev_is_eliminated, int * dev_pivots) {
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		__shared__ int chosen_pivot;
		if(!dev_is_eliminated[v]){
			chosen_pivot = v;
		}

		__syncthreads();

		dev_pivots[v] = chosen_pivot;
	}
}

void update(int num_nodes, int * dev_pivots, bool * fw_is_visited, bool * bw_is_visited, bool * dev_is_eliminated, bool * stop) {
	bool * dev_fw_is_visited, * dev_bw_is_visited, *dev_stop;
	int * dev_colors;
	long * dev_write_id_for_pivots;

	long * write_id_for_pivots = (long*) malloc(4 * num_nodes * sizeof(long));
	int * colors = (int*) malloc(num_nodes * sizeof(int));

	*stop = true;
	for (int i = 0; i < 4 * num_nodes; i++)
        write_id_for_pivots[i] = -1;

	/* for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_UPDATE); */

	HANDLE_ERROR(cudaMalloc((void**)&dev_fw_is_visited, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_bw_is_visited, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_colors, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_write_id_for_pivots, 4 * num_nodes * sizeof(long)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));
	
	
	HANDLE_ERROR(cudaMemcpy(dev_fw_is_visited, fw_is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_bw_is_visited, bw_is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_colors, colors, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_write_id_for_pivots, write_id_for_pivots, 4 * num_nodes * sizeof(long), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_stop, stop, sizeof(bool), cudaMemcpyHostToDevice));
	

	/*
	Dai paper:
	These subgraphs are 
	1) the strongly connected component with the pivot;
	2) the subgraph given by vertices in the forward closure but not in the backward closure; 
	3) the subgraph given by vertices in the backward closure but not in the forward closure;
	4) the subgraph given by vertices that are neither in the forward nor in the backward closure.
	
	The subgraphs that do not contain the pivot form three independent instances of the same problem, and therefore, 
	they are recursively processed in parallel with the same algorithm
	*/
	
	int NOB = num_nodes / 1024 + (num_nodes % 1024 == 0 ? 0 : 1);

	set_colors<<<NOB, 1024>>>(num_nodes, dev_fw_is_visited, dev_bw_is_visited, dev_pivots, dev_colors, dev_is_eliminated, dev_write_id_for_pivots, dev_stop);
	HANDLE_ERROR(cudaMemcpy(fw_is_visited, dev_fw_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(bw_is_visited, dev_bw_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));


	HANDLE_ERROR(cudaFree(dev_fw_is_visited));
	HANDLE_ERROR(cudaFree(dev_bw_is_visited));
	HANDLE_ERROR(cudaFree(dev_stop));
	

	// setto i valori dei pivot che hanno vinto la race	
	set_race_winners<<<NOB, 1024>>>(num_nodes, dev_is_eliminated, dev_pivots, dev_colors, dev_write_id_for_pivots);

	HANDLE_ERROR(cudaMemcpy(colors, dev_colors, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(write_id_for_pivots, dev_write_id_for_pivots, 4 * num_nodes * sizeof(long), cudaMemcpyDeviceToHost));

	HANDLE_ERROR(cudaFree(dev_colors));
	HANDLE_ERROR(cudaFree(dev_write_id_for_pivots));

	/* for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_UPDATE); */

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("colors[" + to_string(i) + "] = ", colors[i], DEBUG_UPDATE);

	for (int i = 0; i < 4 * num_nodes; i++)
        DEBUG_MSG("write_id_for_pivots[" + to_string(i) + "] = ", write_id_for_pivots[i], DEBUG_UPDATE);
}

__global__ void trim_u_kernel(int num_nodes, int * dev_nodes, int * dev_adjacency_list, int * dev_pivots, bool * dev_is_u, int * dev_is_scc){
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if(v < num_nodes){
		if(dev_is_u[v]){
			for(int u = dev_nodes[v]; u < dev_nodes[v+1]; ++u) {
				if(dev_pivots[v] != dev_pivots[dev_adjacency_list[u]]) {
					dev_is_scc[dev_pivots[dev_adjacency_list[u]]] = -1;
				}
			}
		}
	}

}

__global__ void trim_u_propagation(int num_nodes, int * dev_pivots, int * dev_is_scc) {
	int v = threadIdx.x + blockIdx.x * blockDim.x;

	if (v < num_nodes)
		dev_is_scc[v] = dev_is_scc[dev_pivots[v]];
}

__global__ void calculate_more_than_one(int num_nodes, int * more_than_one_dev, int * is_scc_dev) {
	int u = threadIdx.x + blockIdx.x * blockDim.x;

	if (u < num_nodes){
		if(is_scc_dev[u] != -1){
			// printf("%d, more_than_one_dev[%d] = %d\n", u, is_scc_dev[u], more_than_one_dev[is_scc_dev[u]]);
			// atomicAdd può essere migliorato
			atomicAdd(&more_than_one_dev[is_scc_dev[u]], 1);
		}
	}
}

__global__ void is_scc_adjust(int num_nodes, int * more_than_one_dev, int * is_scc_dev) {
	int u = threadIdx.x + blockIdx.x * blockDim.x;

	if (u < num_nodes){
		if(more_than_one_dev[u] == 1)
			is_scc_dev[u] = -1;
	}
}

int count_distinct(int arr[], int n){
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
	// Impostazione del device
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	int num_nodes, num_edges;
    int * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose, * is_scc;
	bool * fw_is_visited, * bw_is_visited, * is_eliminated, * fw_is_expanded, * bw_is_expanded, * is_u;

    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, is_u);

	if(DEBUG_MAIN){
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes[" + to_string(i) + "] = ", nodes[i], DEBUG_MAIN);
		for (int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list[" + to_string(i) + "] = ", adjacency_list[i], DEBUG_MAIN);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("nodes_transpose[" + to_string(i) + "] = ", nodes_transpose[i], DEBUG_MAIN);
		for (int i = 0; i < num_edges; i++)
			DEBUG_MSG("adjacency_list_transpose[" + to_string(i) + "] = ", adjacency_list_transpose[i], DEBUG_MAIN);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("is_u[" + to_string(i) + "] = ", is_u[i], DEBUG_MAIN);
	}

	const int THREADS_PER_BLOCK = prop.maxThreadsPerBlock;
	const int NUMBER_OF_BLOCKS = num_nodes / THREADS_PER_BLOCK + (num_nodes % THREADS_PER_BLOCK == 0 ? 0 : 1);

	// Dichiarazioni di variabili device
	int * dev_is_scc, * dev_more_than_one, * dev_nodes, * dev_adjacency_list, * dev_nodes_transpose, * dev_adjacency_list_transpose, * dev_pivots, *dev_chosen_pivot;
	bool * dev_is_u, * dev_is_eliminated;

	HANDLE_ERROR(cudaMalloc((void**)&dev_nodes, (num_nodes+1) * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_nodes_transpose, (num_nodes+1) * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_adjacency_list, num_edges * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_adjacency_list_transpose, num_edges * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_eliminated, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_pivots, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_chosen_pivot, sizeof(int)));

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpy(dev_nodes, nodes, (num_nodes+1) * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_adjacency_list, adjacency_list, num_edges * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_nodes_transpose, nodes_transpose, (num_nodes+1) * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(int), cudaMemcpyHostToDevice));

	fw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
	bw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
	is_eliminated = (bool*) malloc(num_nodes * sizeof(bool));
	fw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));
	bw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));

	for (int i = 0; i < num_nodes; i++){
		is_eliminated[i] = !is_u[i];
		fw_is_visited[i] = false;
		bw_is_visited[i] = false;
		fw_is_expanded[i] = false;
		bw_is_expanded[i] = false;
	}
	
	HANDLE_ERROR(cudaMemcpy(dev_is_eliminated, is_eliminated, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	trimming(num_nodes, num_edges, dev_nodes, dev_nodes_transpose, dev_adjacency_list, dev_adjacency_list_transpose, dev_is_eliminated, THREADS_PER_BLOCK, NUMBER_OF_BLOCKS);

	initialize_pivot<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_is_eliminated, dev_pivots);
		
    bool stop = false;

    while (!stop){
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, dev_nodes, dev_adjacency_list, dev_pivots, fw_is_visited, dev_is_eliminated, fw_is_expanded);
		
        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, dev_nodes_transpose, dev_adjacency_list_transpose, dev_pivots, bw_is_visited, dev_is_eliminated, bw_is_expanded);

		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, dev_nodes, dev_nodes_transpose, dev_adjacency_list, dev_adjacency_list_transpose, dev_is_eliminated, THREADS_PER_BLOCK, NUMBER_OF_BLOCKS);

		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, dev_pivots, fw_is_visited, bw_is_visited, dev_is_eliminated, &stop);
		
    }
	
    // ---- INIZIO DEBUG ----
	if (DEBUG_FW_BW){
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("fw_is_visited[" + to_string(i) + "] = ", fw_is_visited[i], DEBUG_FW_BW);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("bw_is_visited[" + to_string(i) + "] = ", bw_is_visited[i], DEBUG_FW_BW);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("is_eliminated[" + to_string(i) + "] = ", is_eliminated[i], DEBUG_FW_BW);
		/* for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_FW_BW); */
	}
    // ---- FINE DEBUG ----

	/* for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_MAIN); */

	/*
	Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	Quello che manca è capire quali SCC sono accettabili, ovvero tali che nell'insieme prec(SCC) non ci sia neanche un nodo che appartiene a U
	*/ 

	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1
	// Per iniziare le assegnamo gli stessi valori di pivots, che verranno modificati in seguito
	is_scc = (int*) malloc(num_nodes * sizeof(int));
	
	// Allochiamo more_than_one, che per ogni nodo che fa da pivot viene assegnato un contatore, il quale conta quante volte appare tale pivot
	// Se appare solo una volta, allora il nodo non fa parte di nessuna SCC
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_u, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_scc, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_more_than_one, num_nodes * sizeof(int)));
	
	
	HANDLE_ERROR(cudaMemcpy(dev_is_scc, dev_pivots, num_nodes * sizeof(int), cudaMemcpyDeviceToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_u, is_u, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemset(dev_more_than_one, 0, num_nodes * sizeof(int)));

	trim_u_kernel<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_nodes, dev_adjacency_list, dev_pivots, dev_is_u, dev_is_scc);
	
	trim_u_propagation<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_pivots, dev_is_scc);
	
	calculate_more_than_one<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_more_than_one, dev_is_scc);
	
	is_scc_adjust<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_more_than_one, dev_is_scc);
	
	HANDLE_ERROR(cudaFree(dev_pivots));
	HANDLE_ERROR(cudaFree(dev_is_u));
	HANDLE_ERROR(cudaFree(dev_is_eliminated));
	HANDLE_ERROR(cudaFree(dev_more_than_one));
	HANDLE_ERROR(cudaFree(dev_nodes));
	HANDLE_ERROR(cudaFree(dev_adjacency_list));
	HANDLE_ERROR(cudaFree(dev_nodes_transpose));
	HANDLE_ERROR(cudaFree(dev_adjacency_list_transpose));

	HANDLE_ERROR(cudaMemcpy(is_scc, dev_is_scc, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(dev_is_scc));

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("is_scc[" + to_string(i) + "] = ", is_scc[i], false);


	DEBUG_MSG("Number of SCCs found: ", count_distinct(is_scc, num_nodes), DEBUG_FINAL);
}