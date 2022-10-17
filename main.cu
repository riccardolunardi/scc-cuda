#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>

#include "adj_list.cpp"
#include "reverse_adj_list.cpp"

#define debug 1
using namespace std;

__global__ void trimming(int * nodes, int * nodes_transpose, int * adjacency_list_transpose, int num_nodes, int num_edges, bool * forward_visited, bool * backward_visited, int * subgraph, bool * forward_terminate) {
	int i, list_pointer_1, list_pointer_2;
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	bool elim = true;

	if (id < num_nodes) {
		printf("FW :: v[%d]=%d : v[%d]=%d \n", id, nodes[id], id + 1, nodes[id + 1]);
		//printf("BW :: v[%d]=%d : v[%d]=%d \n", id, nodes_transpose[id], id+1, nodes_transpose[id+1]);

		if (forward_visited[id] == false) {

			if (id == num_nodes - 1) {
				elim = true;
				if (nodes[id] == num_edges || nodes_transpose[id] == num_edges) {
					forward_visited[id] = true;
					backward_visited[id] = true;
					subgraph[id] = 4 * id + 1;
					* forward_terminate = false;
					printf("Trim e : v[%d] : %d, setting forward_terminate to %d \n", id, subgraph[id], * forward_terminate);
				}
			} else if ((nodes[id] == nodes[id + 1]) || (nodes_transpose[id] == nodes_transpose[id + 1])) {
				forward_visited[id] = true;
				backward_visited[id] = true;
				subgraph[id] = 4 * id + 1;
				* forward_terminate = false;
				printf("Trim e : v[%d] : %d, setting forward_terminate to %d \n", id, subgraph[id], * forward_terminate);
			} else {
				list_pointer_1 = nodes_transpose[id];

				if (id == num_nodes - 1)
					list_pointer_2 = num_edges;
				else
					list_pointer_2 = nodes_transpose[id + 1];

				for (i = list_pointer_1; i < list_pointer_2; i++) {
					printf("iteration %d subgraph[%d]=%d, subgraph[%d]=%d\n", i, adjacency_list_transpose[i] - 1, subgraph[adjacency_list_transpose[i] - 1], id, subgraph[id]);
					if (subgraph[adjacency_list_transpose[i] - 1] == subgraph[id]) {
						elim = false;
						break;
					}
				}
				if (elim == true) {
					forward_visited[id] = true;
					backward_visited[id] = true;
					subgraph[id] = 4 * id + 1;
					* forward_terminate = false;
					printf("Trim e : v[%d] : %d, setting forward_terminate to %d,subgraph to %d \n", id, subgraph[id], * forward_terminate, subgraph[id]);
				}

			}
		}
	}
}

__global__ void forward_closure(int * device_nodes, int * device_adjacency_list, int * subgraph, bool * forward_visited, bool * forward_terminate, int num_nodes, int num_edges) {
	//	printf("in fw\n");
	int i, list_pointer_1, list_pointer_2;

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int pivot = id + 1;
	if (id < num_nodes) {
		//		printf("TID = %d v : %d sg %d\n",pivot, forward_visited[id], subgraph[id]);

		if (forward_visited[id]) {
			list_pointer_1 = device_nodes[pivot - 1];

			if (pivot == num_nodes)
				list_pointer_2 = num_edges;
			else
				list_pointer_2 = device_nodes[pivot];

			//			printf("id = %d :: %d %d \n", id, list_pointer_1, list_pointer_2);	
			for (i = list_pointer_1; i < list_pointer_2; i++) {
				//				printf("v[%d] : %d sp=%d s=%d \n", device_adjacency_list[i], forward_visited[device_adjacency_list[i]-1],subgraph[pivot-1],subgraph[device_adjacency_list[i]-1]);	
				if (forward_visited[device_adjacency_list[i] - 1] == false && subgraph[pivot - 1] == subgraph[device_adjacency_list[i] - 1]) {
					//					printf("src -> dest : %d -> %d\n",pivot, device_adjacency_list[i]);
					forward_visited[device_adjacency_list[i] - 1] = true;
					* forward_terminate = false;
				}
			}
		}
	}
}

__global__ void generate_subgraph(int pivot, bool * forward_visited, bool * backward_visited, int * subgraph, int num_nodes) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < num_nodes) {
		// Il nodo è stato visitato sia dalla backward che dalla forward
		// Il nodo fa parte di una SCC
		if (forward_visited[id] == backward_visited[id] && forward_visited[id] == true) {
			subgraph[id] = 4 * pivot;
		}

		// Il nodo "id" è stato visitato dalla forward, ma non dalla backward
		// Si deve calcolare FB(F\B)
		if (forward_visited[id] != backward_visited[id] && forward_visited[id] == true) {
			subgraph[id] = 4 * pivot + 1;
			forward_visited[id] = backward_visited[id] = false;
		}

		// Il nodo "id" è stato visitato dalla backward, ma non dalla forward
		// Si deve calcolare FB(B\F)
		if (forward_visited[id] != backward_visited[id] && backward_visited[id] == true) {
			subgraph[id] = 4 * pivot + 2;
			forward_visited[id] = backward_visited[id] = false;
		}

		// Il nodo non è stato visitato da nessuno
		// Si deve calcolare FB( V \ (B U F))
		if (forward_visited[id] == backward_visited[id] && backward_visited[id] == false) {
			subgraph[id] = 4 * pivot + 3;
			forward_visited[id] = backward_visited[id] = false;
		}
	}
}

void fw_bw(int num_nodes, int num_edges, int * nodes, int * adjacency_list, int * nodes_transpose, int * adjacency_list_transpose) {
	int * device_nodes, * device_adjacency_list, * device_nodes_transpose, * device_adjacency_list_transpose;
	int * subgraph;
	bool * forward_visited, * backward_visited;
	bool * forward_terminate, * backward_terminate, * device_forward_terminate, * device_backward_terminate;
	int i = 0;

	cudaMalloc((void ** ) & device_nodes, num_nodes * (sizeof(int)));
	cudaMalloc((void ** ) & device_adjacency_list, num_edges * (sizeof(int)));
	cudaMalloc((void ** ) & device_nodes_transpose, num_nodes * (sizeof(int)));
	cudaMalloc((void ** ) & device_adjacency_list_transpose, num_edges * (sizeof(int)));
	cudaMalloc((void ** ) & subgraph, num_nodes * (sizeof(int)));
	cudaMalloc((void ** ) & forward_visited, num_nodes * (sizeof(bool)));
	cudaMalloc((void ** ) & backward_visited, num_nodes * (sizeof(bool)));

	cudaHostAlloc((void ** ) & forward_terminate, 1 * sizeof(bool), cudaHostAllocMapped);
	cudaHostAlloc((void ** ) & backward_terminate, 1 * sizeof(bool), cudaHostAllocMapped);

	cudaMemset(subgraph, 0, num_nodes * sizeof(int));
	cudaMemset(forward_visited, false, num_nodes);
	cudaMemset(backward_visited, false, num_nodes);

	cudaHostGetDevicePointer( & device_forward_terminate, forward_terminate, 0);
	cudaHostGetDevicePointer( & device_backward_terminate, backward_terminate, 0);

	cudaMemcpy(device_nodes, nodes, sizeof(int) * num_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(device_adjacency_list, adjacency_list, sizeof(int) * num_edges, cudaMemcpyHostToDevice);
	cudaMemcpy(device_nodes_transpose, nodes_transpose, sizeof(int) * num_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(device_adjacency_list_transpose, adjacency_list_transpose, sizeof(int) * num_edges, cudaMemcpyHostToDevice);

	int num_blocks, num_threads_per_block, pivot;
	num_threads_per_block = 256;
	num_blocks = num_nodes / num_threads_per_block + (num_nodes % num_threads_per_block == 0 ? 0 : 1);

	if (debug)
		cout << "N° blocks: " << num_blocks << " ,n° threads: " << num_threads_per_block << endl;

	// Complete Trimming
	// while (*forward_terminate == false || i < 5) {                  // Il *forward_terminate == false era commentato, ho commentato tutta la riga invece
    // DOMANDA: perchè fa 'trimming 5 volte?'
	while (i < 5) {
		* forward_terminate = true;
		trimming << < num_blocks, num_threads_per_block >>> (device_nodes, device_nodes_transpose, device_adjacency_list_transpose, num_nodes, num_edges, forward_visited, backward_visited, subgraph, device_forward_terminate);
		cudaThreadSynchronize();
		printf("terminate : %d \n", * forward_terminate);
		i++;
	}

	* forward_terminate = false;
	pivot = 0;
	cudaMemset( & forward_visited[pivot], true, 1);
	cudaMemset( & backward_visited[pivot], true, 1);

	//Forward-Closure
	if (debug) cout << "Forward closure\n";
	while ( * forward_terminate == false) {
		* forward_terminate = true;
		forward_closure << < num_blocks, num_threads_per_block >>> (device_nodes, device_adjacency_list, subgraph, forward_visited, device_forward_terminate, num_nodes, num_edges);
		cudaThreadSynchronize();
	}

	//Backward-Closure
	if (debug) cout << "Backward  closure\n";
	while ( * backward_terminate == false) {
		* backward_terminate = true;
		forward_closure << < num_blocks, num_threads_per_block >>> (device_nodes_transpose, device_adjacency_list_transpose, subgraph, backward_visited, device_backward_terminate, num_nodes, num_edges);
		cudaThreadSynchronize();
	}

	//Finding 4 Subgraphs		
	generate_subgraph << < num_blocks, num_threads_per_block >>> (pivot, forward_visited, backward_visited, subgraph, num_nodes);

	cudaFree(device_nodes);
	cudaFree(device_adjacency_list);
	cudaFree(device_nodes_transpose);
	cudaFree(device_adjacency_list_transpose);
	cudaFree(subgraph);
	cudaFree(forward_visited);
	cudaFree(backward_visited);
	cudaFreeHost(forward_terminate);
	cudaFreeHost(backward_terminate);
}

int main(int argc, char ** argv) {
	if (argc != 2) {
		cout << " Invalid Usage !! Usage is ./a.out <graph_input_file> \n";
		return -1;
	}
	const char *filename = argv[1];

	//Così sembra che la prima riga sia letteralmente scartata
	std::string line;
	std::ifstream infile(filename);
	std::getline(infile, line);

	/*
	La seconda riga contiene:
	- Numero di archi
	- Numero di nodi
	- Numero bi BOH
	*/
	std::getline(infile, line);
	std::istringstream iss(line);

	char percentage_sign;
	int num_edges, num_nodes, edge_weight, i;

	iss >> percentage_sign; //Questo dovrebbe essere il simbolo "%" nei file, che viene acquisito qui, ma poi mai più usato
	iss >> num_edges;
	iss >> num_nodes;
	iss >> edge_weight;

	infile.close();

	std::cout << "Number of nodes: " << num_nodes << endl;
	std::cout << "Number of edges: " << num_edges << endl;

	//Inizializzazione delle strutture dati principali
	int *nodes = new int[num_nodes];
	int *adjacency_list = new int[num_edges];
	int *nodes_transpose = new int[num_nodes];
	int *adjacency_list_transpose = new int[num_edges];

    // Inizializzazione delle liste
	for (i = 0; i < num_nodes; i++){
		nodes[i] = 0;
		nodes_transpose[i] = 0;
	}

	// Creazione delle liste di adiacenza
	adj_list(filename, nodes, adjacency_list);

	// Creazione delle liste di adiacenza del grafo trasposto (per la backward clousure)
	// Forse si può evitare la ripetizione di codice usando il codice di adj_list leggermente modificato
	reverse_adj_list(filename, nodes_transpose, adjacency_list_transpose);

	if (debug) {
		cout << " Adj List " << endl;
		cout << " ---O(V) \n";
		for (i = 0; i < num_nodes; i++) {
			cout << "nodes[" << i << "] : " << nodes[i] << endl;
		}
		cout << " ---O(E) \n";
		for (i = 0; i < num_edges; i++) {
			cout << "adjacency_list[" << i << "] : " << adjacency_list[i] << endl;
		}

		cout << " Transpose Adj List " << endl;
		cout << " ---O(V) \n";
		for (i = 0; i < num_nodes; i++) {
			cout << "nodes[" << i << "] : " << nodes_transpose[i] << endl;
		}
		cout << " ---O(E) \n";
		for (i = 0; i < num_edges; i++) {
			cout << "adjacency_list[" << i << "] : " << adjacency_list_transpose[i] << endl;
		}
	}

	fw_bw(num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose);

	delete(nodes);
	delete(adjacency_list);
	delete(nodes_transpose);
	delete(adjacency_list_transpose);
}