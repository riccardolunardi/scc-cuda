#include "utils.cpp"
#include <cstring>
#include <cuda.h>
using namespace std;

#define DEBUG_F_KERNEL false
#define DEBUG_REACH true
#define DEBUG_TRIMMING_KERNEL false
#define DEBUG_TRIMMING true
#define DEBUG_UPDATE true
#define DEBUG_FW_BW true
#define DEBUG_MAIN true
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
	//printf("Checking %d\n", v);
    // per ogni nodo
	if(v < num_nodes) {
        //printf("dev_is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
		//printf("dev_is_visited[%d] -> %d\n", v, dev_is_visited[v]);
		//printf("dev_is_expanded[%d] -> %d\n", v, dev_is_expanded[v]);
        /* DEBUG_MSG("is_visited[" << v << "] -> ", dev_is_visited[v], DEBUG_F_KERNEL);
        DEBUG_MSG("is_expanded["<< v <<"] -> ", dev_is_expanded[v], DEBUG_F_KERNEL); */
		
        // si controlla se non è stato eliminato E è stato eliminato E non è stato espanso
		if(!dev_is_eliminated[v] && dev_is_visited[v] && !dev_is_expanded[v]) {
            // si segna come espanso
			dev_is_expanded[v] = true;
			//printf("	u va da %d a %d\n", dev_nodes[v], dev_nodes[v+1]);

            // per ogni nodo a cui punta
			int row_begin = dev_nodes[v];
			int row_end = dev_nodes[v + 1];
			for(int u = row_begin; u < row_end; u++) {	
					
				/* printf("		Nodo %d connesso a nodo %d\n", v, dev_adjacency_list[u]);	
				printf("		is_eliminated[%d] -> %d\n", dev_adjacency_list[u], dev_is_eliminated[dev_adjacency_list[u]]);
				printf("		is_visited[%d] -> %d\n", dev_adjacency_list[u], dev_is_visited[dev_adjacency_list[u]]);
				printf("		pivots[%d] == pivots[%d] -> %d == %d\n", v, dev_adjacency_list[u],  dev_pivots[v], dev_pivots[dev_adjacency_list[u]]); */

                // si controlla se non è stato eliminato E se non è stato visitato E se il colore del nodo che punta corrisponde a quello del nodo puntato
				int dst = dev_adjacency_list[u];
				
				if(!dev_is_eliminated[dst] && !dev_is_visited[dst] && dev_pivots[v] == dev_pivots[dst]) {
					// printf("			is_visited[%d] -> TRUE\n", dev_adjacency_list[u]);
                    // setta il nodo puntato a visitato
					dev_is_visited[dst] = true;
                    // permette di continuare il ciclo in reach, perchè si è trovato un altro nodo da visitare
					*dev_stop = false;
					//printf("			%d\n", *dev_stop);
				}
			}
		}
	}
}

void reach(int num_nodes, int num_edges, int * dev_nodes, int * dev_adjacency_list, int * pivots, bool * is_visited, bool * is_eliminated, bool * is_expanded) {
    // Tutti i pivot vengono segnati come visitati
    //printf("Ciao %d\n", 4);

	for(int i=0; i < num_nodes; i++) {
        is_visited[ pivots[i] ] = true;
    }

	//printf("Ciao %d\n", 5);

	int *dev_pivots;
	bool *dev_is_visited, *dev_is_eliminated, *dev_is_expanded;
	bool stop, *dev_stop;
	stop = false;
	//printf("Ciao %d\n", 5);

	int NOB = num_nodes / 1024 + (num_nodes % 1024 == 0 ? 0 : 1);

	HANDLE_ERROR(cudaMalloc((void**)&dev_pivots, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_visited, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_eliminated, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_expanded, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));
	cudaDeviceSynchronize();

	HANDLE_ERROR(cudaMemcpy(dev_pivots, pivots, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_visited, is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_eliminated, is_eliminated, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_expanded, is_expanded, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();

    // si effettua la chiusura in avanti/indietro
    while(!stop) {
		stop = true;
		HANDLE_ERROR(cudaMemcpy(dev_stop, &stop, sizeof(bool), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();
		
        f_kernel<<<NOB, 1024>>>(num_nodes, num_edges, dev_nodes, dev_adjacency_list, dev_pivots, dev_is_visited, dev_is_eliminated, dev_is_expanded, dev_stop);
		cudaDeviceSynchronize();
		for (int i = 0; i < num_nodes; i++){
			DEBUG_MSG("is_visited[" << i << "] -> ", is_visited[i], DEBUG_REACH);
		}

		HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
		cudaDeviceSynchronize();
		//printf("stop: %d\n", stop);
    }

	HANDLE_ERROR(cudaMemcpy(is_visited, dev_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(is_eliminated, dev_is_eliminated, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(is_expanded, dev_is_expanded, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(pivots, dev_pivots, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();


	HANDLE_ERROR(cudaFree(dev_pivots));
	HANDLE_ERROR(cudaFree(dev_is_visited));
	HANDLE_ERROR(cudaFree(dev_is_eliminated));
	HANDLE_ERROR(cudaFree(dev_is_expanded));
	HANDLE_ERROR(cudaFree(dev_stop));
	cudaDeviceSynchronize();

}

__global__ void trimming_kernel(int num_nodes, int num_edges, int * dev_nodes, int * dev_nodes_transpose, int * dev_adjacency_list, int * dev_pivots, bool * dev_is_eliminated, bool * dev_stop){
	int v = threadIdx.x + blockIdx.x * blockDim.x;
	bool elim;

	if(v < num_nodes) {
		//printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
		if(!dev_is_eliminated[v]){
			elim = true;
			
			// Nel caso un nodo abbia entrambi in_degree o out_degree diversi da 0 allora non va eliminato
			if(dev_nodes[v] != dev_nodes[v+1] && dev_nodes_transpose[v] != dev_nodes_transpose[v+1]){
				elim = false;
			}

			// Nel caso un arco di v faccia parte dello stesso sottografo, allora non va eliminato
			// Non serve farlo anche per la lista trasposta perchè alla fine l'if sui pivots e sarebbe la stessa cosa
			//printf("	nodo %d e nodo %d\n", v, v+1);
			//printf("	u va da %d a %d\n", dev_nodes[v], dev_nodes[v+1]);
			for(int u = dev_nodes[v]; u < dev_nodes[v+1]; u++){
			 	//printf("adjacency_list[%d] -> %d\n", u, dev_adjacency_list[u]);
			 	// elim potrebbe ovviamente essere messo prima del for, per evitare cicli in più
				// forse va mosso, in base a come capiremo ottimizzare in CUDA
				if(dev_pivots[dev_adjacency_list[u]] == dev_pivots[v] && !elim){
					//printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
			 		elim = false;
				}	
			}

			if(elim){
				dev_is_eliminated[v] = true;
				//printf("is_eliminated[%d] -> %d\n", v, dev_is_eliminated[v]);
				*dev_stop = false;
			}
		}
	}
}

void trimming(int num_nodes, int num_edges, int * dev_nodes, int * dev_nodes_transpose, int * dev_adjacency_list, int * pivots, bool * is_eliminated) {
	int * dev_pivots;
	bool * dev_is_eliminated, * dev_stop, stop;
	stop = false;

	HANDLE_ERROR(cudaMalloc((void**)&dev_pivots, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_eliminated, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));
	int NOB = num_nodes / 1024 + (num_nodes % 1024 == 0 ? 0 : 1);

	HANDLE_ERROR(cudaMemcpy(dev_pivots, pivots, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_eliminated, is_eliminated, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_stop, &stop, sizeof(bool), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	}
    while(!stop) {
        stop = true;
		cudaDeviceSynchronize();
		HANDLE_ERROR(cudaMemcpy(dev_stop, &stop, sizeof(bool), cudaMemcpyHostToDevice));
		cudaDeviceSynchronize();
        trimming_kernel<<<NOB, 1024>>>(num_nodes, num_edges, dev_nodes, dev_nodes_transpose, dev_adjacency_list, dev_pivots, dev_is_eliminated, dev_stop);
		cudaDeviceSynchronize();
		

		HANDLE_ERROR(cudaMemcpy(&stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
		//printf("stop: %d\n", stop);
		cudaDeviceSynchronize();
    }

	HANDLE_ERROR(cudaMemcpy(is_eliminated, dev_is_eliminated, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();
	HANDLE_ERROR(cudaFree(dev_is_eliminated));
	HANDLE_ERROR(cudaFree(dev_pivots));
	HANDLE_ERROR(cudaFree(dev_stop));
	cudaDeviceSynchronize();
	for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("is_eliminated[" << i << "] -> ", is_eliminated[i], DEBUG_TRIMMING);
	}
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

void update(int num_nodes, int * pivots, bool * fw_is_visited, bool * bw_is_visited, bool * is_eliminated, bool * stop) {
	bool * dev_fw_is_visited, * dev_bw_is_visited, * dev_is_eliminated, *dev_stop;
	int * dev_colors, *dev_pivots;
	long * dev_write_id_for_pivots;

	long * write_id_for_pivots = (long*) malloc(4 * num_nodes * sizeof(long));
	int * colors = (int*) malloc(num_nodes * sizeof(int));

	*stop = true;
	for (int i = 0; i < 4 * num_nodes; i++)
        write_id_for_pivots[i] = -1;

	HANDLE_ERROR(cudaMalloc((void**)&dev_fw_is_visited, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_bw_is_visited, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_is_eliminated, num_nodes * sizeof(bool)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_colors, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_pivots, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_write_id_for_pivots, 4 * num_nodes * sizeof(long)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_stop, sizeof(bool)));
	cudaDeviceSynchronize();
	
	HANDLE_ERROR(cudaMemcpy(dev_fw_is_visited, fw_is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_bw_is_visited, bw_is_visited, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_is_eliminated, is_eliminated, num_nodes * sizeof(bool), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_colors, colors, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_pivots, pivots, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_write_id_for_pivots, write_id_for_pivots, 4 * num_nodes * sizeof(long), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_stop, stop, sizeof(bool), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();

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
	cudaDeviceSynchronize();
	HANDLE_ERROR(cudaMemcpy(fw_is_visited, dev_fw_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(bw_is_visited, dev_bw_is_visited, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(is_eliminated, dev_is_eliminated, num_nodes * sizeof(bool), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(colors, dev_colors, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(pivots, dev_pivots, num_nodes * sizeof(int), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(write_id_for_pivots, dev_write_id_for_pivots, 4 * num_nodes * sizeof(long), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(stop, dev_stop, sizeof(bool), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();
	//printf("STOP? %d\n", *stop);

	HANDLE_ERROR(cudaFree(dev_fw_is_visited));
	HANDLE_ERROR(cudaFree(dev_bw_is_visited));
	HANDLE_ERROR(cudaFree(dev_is_eliminated));
	HANDLE_ERROR(cudaFree(dev_colors));
	HANDLE_ERROR(cudaFree(dev_pivots));
	HANDLE_ERROR(cudaFree(dev_write_id_for_pivots));
	HANDLE_ERROR(cudaFree(dev_stop));
	cudaDeviceSynchronize();
	
	for (int i = 0; i < 4 * num_nodes; i++)
        DEBUG_MSG("write_id_for_pivots[" + to_string(i) + "] = ", write_id_for_pivots[i], true);

	// setto i valori dei pivot che hanno vinto la race
	// in CUDA questo è da fare dopo una sincronizzazione
	cudaDeviceSynchronize();
	for (int i = 0; i < num_nodes; i++) {
		if(is_eliminated[i]){
			pivots[i] = i;
		}else{
			pivots[i] = write_id_for_pivots[colors[i]];
		}
	}

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_UPDATE);
}

void trim_u_kernel(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * pivots, bool * is_u, int *& is_scc) {
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

void trim_u_propagation_h(int num_nodes, int * pivots, int *& is_scc) {
	for(int u = 0; u < num_nodes; ++u ) {
		is_scc[u] = is_scc[pivots[u]];
	}
}

__global__ void trim_u_propagation(int num_nodes, int * pivots, int * is_scc) {
	int u = threadIdx.x + blockIdx.x * blockDim.x;

	if (u < num_nodes)
		is_scc[u] = is_scc[pivots[u]];
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
	// Impostazione del device e dei flag
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);
	cudaSetDeviceFlags(cudaDeviceMapHost);

    if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> \n";
		return -1;
	}

	int num_nodes, num_edges;
    int * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose, * pivots, * is_scc, * more_than_one;
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
	int * dev_is_scc, * dev_more_then_one, * dev_nodes, * dev_adjacency_list, * dev_nodes_transpose, * dev_adjacency_list_transpose;

	HANDLE_ERROR(cudaMalloc((void**)&dev_nodes, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_nodes_transpose, num_nodes * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_adjacency_list, num_edges * sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void**)&dev_adjacency_list_transpose, num_edges * sizeof(int)));

	fw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
	bw_is_visited = (bool*) malloc(num_nodes * sizeof(bool));
	is_eliminated = (bool*) malloc(num_nodes * sizeof(bool));
	fw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));
	bw_is_expanded = (bool*) malloc(num_nodes * sizeof(bool));
	pivots = (int*) malloc(num_nodes * sizeof(int));

	for (int i = 0; i < num_nodes; i++){
		is_eliminated[i] = !is_u[i];
		pivots[i] = 5;
	}

	memset(fw_is_visited, false, num_nodes * sizeof(bool));
	memset(bw_is_visited, false, num_nodes * sizeof(bool));
	memset(fw_is_expanded, false, num_nodes * sizeof(bool));
	memset(bw_is_expanded, false, num_nodes * sizeof(bool));

	/* for (int i = 0; i < num_nodes; i++){
		DEBUG_MSG("fw_is_visited[" + to_string(i) + "] = ", fw_is_visited[i], DEBUG_MAIN);
		DEBUG_MSG("bw_is_visited[" + to_string(i) + "] = ", bw_is_visited[i], DEBUG_MAIN);
		DEBUG_MSG("fw_is_expanded[" + to_string(i) + "] = ", fw_is_expanded[i], DEBUG_MAIN);
		DEBUG_MSG("bw_is_expanded[" + to_string(i) + "] = ", bw_is_expanded[i], DEBUG_MAIN);
		DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_MAIN);
	} */

	// Le strutture principali le copiamo nel device già qui, visto che non verranno mai modificate
	HANDLE_ERROR(cudaMemcpy(dev_nodes, nodes, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_adjacency_list, adjacency_list, num_edges * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_nodes_transpose, nodes_transpose, num_nodes * sizeof(int), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(dev_adjacency_list_transpose, adjacency_list_transpose, num_edges * sizeof(int), cudaMemcpyHostToDevice));

    bool stop = false;

    while (!stop){
		DEBUG_MSG("Forward reach:" , "", DEBUG_FW_BW);
        reach(num_nodes, num_edges, dev_nodes, dev_adjacency_list, pivots, fw_is_visited, is_eliminated, fw_is_expanded);
		cudaDeviceSynchronize();

        DEBUG_MSG("Backward reach:" , "", DEBUG_FW_BW);
		reach(num_nodes, num_edges, dev_nodes_transpose, dev_adjacency_list_transpose, pivots, bw_is_visited, is_eliminated, bw_is_expanded);
		cudaDeviceSynchronize();

		DEBUG_MSG("Trimming:" , "", DEBUG_FW_BW);
        trimming(num_nodes, num_edges, dev_nodes, dev_nodes_transpose, dev_adjacency_list_transpose, pivots, is_eliminated);
		cudaDeviceSynchronize();

		DEBUG_MSG("Update:" , "", DEBUG_FW_BW);
		update(num_nodes, pivots, fw_is_visited, bw_is_visited, is_eliminated, &stop);
		cudaDeviceSynchronize();
    }
	cudaDeviceSynchronize();
    // ---- INIZIO DEBUG ----
	if (DEBUG_FW_BW){
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("fw_is_visited[" + to_string(i) + "] = ", fw_is_visited[i], DEBUG_FW_BW);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("bw_is_visited[" + to_string(i) + "] = ", bw_is_visited[i], DEBUG_FW_BW);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("is_eliminated[" + to_string(i) + "] = ", is_eliminated[i], DEBUG_FW_BW);
		for (int i = 0; i < num_nodes; i++)
			DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_FW_BW);
	}
    // ---- FINE DEBUG ----

	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("pivots[" + to_string(i) + "] = ", pivots[i], DEBUG_MAIN);

	/*
	Tramite fw_bw_ abbiamo ottenuto, per ogni nodo, il pivot della SCC a cui appartiene.
	Quello che manca è capire quali SCC sono accettabili, ovvero tali che nell'insieme prec(SCC) non ci sia neanche un nodo che appartiene a U
	*/ 


	// Allochiamo is_scc, che alla fine avrà per ogni nodo il pivot della sua SCC se la sua SCC è accettabile, altrimenti -1
	// Per iniziare le assegnamo gli stessi valori di pivots, che verranno modificati in seguito
	HANDLE_ERROR(cudaHostAlloc((void**)&is_scc, num_nodes * sizeof(int), cudaHostAllocMapped));
	HANDLE_ERROR(cudaMemcpy(is_scc, pivots, num_nodes * sizeof(int), cudaMemcpyHostToHost));

	// TODO: Parallelizzare queste due funzioni
	trim_u_kernel(num_nodes, num_edges, nodes, adjacency_list, pivots, is_u, is_scc);
	trim_u_propagation_h(num_nodes, pivots, is_scc);
	

	// Allochiamo more_than_one, che per ogni nodo che fa da pivot viene assegnato un contatore, il quale conta quante volte appare tale pivot
	// Se appare solo una volta, allora il nodo non fa parte di nessuna SCC
	HANDLE_ERROR(cudaHostAlloc((void**)&more_than_one, num_nodes * sizeof(int), cudaHostAllocMapped));
	HANDLE_ERROR(cudaMemset(more_than_one, 0, num_nodes * sizeof(int)));
	
	cudaHostGetDevicePointer(&dev_is_scc, is_scc, 0);
	cudaHostGetDevicePointer(&dev_more_then_one, more_than_one, 0);

	calculate_more_than_one<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_more_then_one, dev_is_scc);
	cudaDeviceSynchronize();
	is_scc_adjust<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(num_nodes, dev_more_then_one, dev_is_scc);
	cudaDeviceSynchronize();
	
	for (int i = 0; i < num_nodes; i++)
        DEBUG_MSG("is_scc[" + to_string(i) + "] = ", is_scc[i], false);

	DEBUG_MSG("Number of SCCs found: ", count_distinct(is_scc, num_nodes), true);
}