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

// cioè almeno sembra, ma se divide in 4 sta facendo OBF? //cristo sandro

__global__ void trimming(int * Vertices, int * VerticesT, int * AdjListT, int n, int m, bool * VisitedF, bool * VisitedB, int * subgraph, bool * terminateF) {
	int i, ListPointer1, ListPointer2;
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	bool elim = true;

	if (id < n) {
		printf("FW :: v[%d]=%d : v[%d]=%d \n", id, Vertices[id], id + 1, Vertices[id + 1]);
		//printf("BW :: v[%d]=%d : v[%d]=%d \n", id, VerticesT[id], id+1, VerticesT[id+1]);

		if (VisitedF[id] == false) {

			if (id == n - 1) {
				elim = true;
				if (Vertices[id] == m || VerticesT[id] == m) {
					VisitedF[id] = true;
					VisitedB[id] = true;
					subgraph[id] = 4 * id + 1;
					* terminateF = false;
					printf("Trim e : v[%d] : %d, setting terminateF to %d \n", id, subgraph[id], * terminateF);
				}
			} else if ((Vertices[id] == Vertices[id + 1]) || (VerticesT[id] == VerticesT[id + 1])) {
				VisitedF[id] = true;
				VisitedB[id] = true;
				subgraph[id] = 4 * id + 1;
				* terminateF = false;
				printf("Trim e : v[%d] : %d, setting terminateF to %d \n", id, subgraph[id], * terminateF);
			} else {
				ListPointer1 = VerticesT[id];

				if (id == n - 1)
					ListPointer2 = m;
				else
					ListPointer2 = VerticesT[id + 1];

				for (i = ListPointer1; i < ListPointer2; i++) {
					printf("iteration %d subgraph[%d]=%d, subgraph[%d]=%d\n", i, AdjListT[i] - 1, subgraph[AdjListT[i] - 1], id, subgraph[id]);
					if (subgraph[AdjListT[i] - 1] == subgraph[id]) {
						elim = false;
						break;
					}
				}
				if (elim == true) {
					VisitedF[id] = true;
					VisitedB[id] = true;
					subgraph[id] = 4 * id + 1;
					* terminateF = false;
					printf("Trim e : v[%d] : %d, setting terminateF to %d,subgraph to %d \n", id, subgraph[id], * terminateF, subgraph[id]);
				}

			}
		}
	}
}

__global__ void forward_closure(int * dVertices, int * dAdjList, int * subgraph, bool * visitedF, bool * terminateF, int numVertices, int numEdges) {
	//	printf("in fw\n");
	int i, ListPointer1, ListPointer2;

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	int pivot = id + 1;
	if (id < numVertices) {
		//		printf("TID = %d v : %d sg %d\n",pivot, visitedF[id], subgraph[id]);

		if (visitedF[id]) {
			ListPointer1 = dVertices[pivot - 1];

			if (pivot == numVertices)
				ListPointer2 = numEdges;
			else
				ListPointer2 = dVertices[pivot];

			//			printf("id = %d :: %d %d \n", id, ListPointer1, ListPointer2);	
			for (i = ListPointer1; i < ListPointer2; i++) {
				//				printf("v[%d] : %d sp=%d s=%d \n", dAdjList[i], visitedF[dAdjList[i]-1],subgraph[pivot-1],subgraph[dAdjList[i]-1]);	
				if (visitedF[dAdjList[i] - 1] == false && subgraph[pivot - 1] == subgraph[dAdjList[i] - 1]) {
					//					printf("src -> dest : %d -> %d\n",pivot, dAdjList[i]);
					visitedF[dAdjList[i] - 1] = true;
					* terminateF = false;
				}
			}
		}
	}
}

__global__ void generate_subgraph(int pivot, bool * visitedF, bool * visitedB, int * subgraph, int numVertices) {
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < numVertices) {
		// Il nodo è stato visitato sia dalla backward che dalla forward
		// Il nodo fa parte di una SCC
		if (visitedF[id] == visitedB[id] && visitedF[id] == true) {
			subgraph[id] = 4 * pivot;
		}

		// Il nodo "id" è stato visitato dalla forward, ma non dalla backward
		// Si deve calcolare FB(F\B)
		if (visitedF[id] != visitedB[id] && visitedF[id] == true) {
			subgraph[id] = 4 * pivot + 1;
			visitedF[id] = visitedB[id] = false;
		}

		// Il nodo "id" è stato visitato dalla backward, ma non dalla forward
		// Si deve calcolare FB(B\F)
		if (visitedF[id] != visitedB[id] && visitedB[id] == true) {
			subgraph[id] = 4 * pivot + 2;
			visitedF[id] = visitedB[id] = false;
		}

		// Il nodo non è stato visitato da nessuno
		// Si deve calcolare FB( V \ (B U F))
		if (visitedF[id] == visitedB[id] && visitedB[id] == false) {
			subgraph[id] = 4 * pivot + 3;
			visitedF[id] = visitedB[id] = false;
		}
	}
}

void fw_bw(int n, int m, int * Vertices, int * AdjacencyList, int * Vertices_Transpose, int * AdjacencyList_Transpose) {
	int * dVertices, * dAdjList, * dVerticesT, * dAdjListT;
	int * subgraph;
	bool * visitedF, * visitedB;
	bool * terminateF, * terminateB, * dterminateF, * dterminateB;
	int i = 0;

	cudaMalloc((void ** ) & dVertices, n * (sizeof(int)));
	cudaMalloc((void ** ) & dAdjList, m * (sizeof(int)));
	cudaMalloc((void ** ) & dVerticesT, n * (sizeof(int)));
	cudaMalloc((void ** ) & dAdjListT, m * (sizeof(int)));
	cudaMalloc((void ** ) & subgraph, n * (sizeof(int)));
	cudaMalloc((void ** ) & visitedF, n * (sizeof(bool)));
	cudaMalloc((void ** ) & visitedB, n * (sizeof(bool)));

	cudaHostAlloc((void ** ) & terminateF, 1 * sizeof(bool), cudaHostAllocMapped);
	cudaHostAlloc((void ** ) & terminateB, 1 * sizeof(bool), cudaHostAllocMapped);

	cudaMemset(subgraph, 0, n * sizeof(int));
	cudaMemset(visitedF, false, n);
	cudaMemset(visitedB, false, n);

	cudaHostGetDevicePointer( & dterminateF, terminateF, 0);
	cudaHostGetDevicePointer( & dterminateB, terminateB, 0);

	cudaMemcpy(dVertices, Vertices, sizeof(int) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(dAdjList, AdjacencyList, sizeof(int) * m, cudaMemcpyHostToDevice);
	cudaMemcpy(dVerticesT, Vertices_Transpose, sizeof(int) * n, cudaMemcpyHostToDevice);
	cudaMemcpy(dAdjListT, AdjacencyList_Transpose, sizeof(int) * m, cudaMemcpyHostToDevice);

	int numBlocks, numThreadsPerBlock, pivot;
	numThreadsPerBlock = 256;
	numBlocks = n / numThreadsPerBlock + (n % numThreadsPerBlock == 0 ? 0 : 1);

	if (debug)
		cout << "N° blocks: " << numBlocks << " ,n° threads: " << numThreadsPerBlock << endl;

	// Complete Trimming
	FWD reach c'è (anche bwd), trimming c'è pivot non ho visto come sia
	// Quindi prima di tutto chiama trimming, yess
	while ( /**terminateF == false ||*/ i < 5) {
		* terminateF = true;
		trimming << < numBlocks, numThreadsPerBlock >>> (dVertices, dVerticesT, dAdjListT, n, m, visitedF, visitedB, subgraph, dterminateF);
		cudaThreadSynchronize();
		printf("terminate : %d \n", * terminateF);
		i++;
	}

	* terminateF = false;
	pivot = 0; //bella
	cudaMemset( & visitedF[pivot], true, 1);
	cudaMemset( & visitedB[pivot], true, 1);

	//Forward-Closure
	if (debug) cout << "Forward closure\n";
	while ( * terminateF == false) {
		* terminateF = true;
		forward_closure << < numBlocks, numThreadsPerBlock >>> (dVertices, dAdjList, subgraph, visitedF, dterminateF, n, m);
		cudaThreadSynchronize();
	}

	//Backward-Closure
	if (debug) cout << "Backward  closure\n";
	while ( * terminateB == false) {
		* terminateB = true;
		forward_closure << < numBlocks, numThreadsPerBlock >>> (dVerticesT, dAdjListT, subgraph, visitedB, dterminateB, n, m);
		cudaThreadSynchronize();
	}

	//Finding 4 Subgraphs		
	generate_subgraph << < numBlocks, numThreadsPerBlock >>> (pivot, visitedF, visitedB, subgraph, n);

	// Tu ti ricordi come si fa la ricorsione in cuda? io no
	// forse va fatto con la programmazione dinamica?
	cudaFree(dVertices);
	cudaFree(dAdjList);
	cudaFree(dVerticesT);
	cudaFree(dAdjListT);
	cudaFree(subgraph);
	cudaFree(visitedF);
	cudaFree(visitedB);
	cudaFreeHost(terminateF);
	cudaFreeHost(terminateB);
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

	char c;
	int m, n, x, i;

	iss >> c; //Questo dovrebbe essere il simbolo "%" nei file, che viene acquisito qui, ma poi mai più usato
	iss >> m;
	iss >> n;
	iss >> x;

	std::cout << "Number of vertices: " << n << endl;
	std::cout << "Number of edges: " << m << endl;

	//Inizializzazione delle strutture dati principali
	int *Vertices = new int[n];
	int *AdjacencyList = new int[m];

	int *Vertices_Transpose = new int[n];
	int *AdjacencyList_Transpose = new int[m];

	for (i = 0; i < n; i++){
		Vertices[i] = 0;
		Vertices_Transpose[i] = 0;
	}

	infile.close();

	// Creazione delle liste di adiacenza
	adj_list(filename, Vertices, AdjacencyList);

	// Creazione delle liste di adiacenza del grafo trasposto (per la backward clousure)
	// Forse si può evitare la ripetizione di codice usando il codice di adj_list leggermente modificato
	reverse_adj_list(filename, Vertices_Transpose, AdjacencyList_Transpose);

	if (debug) {
		cout << " Adj List " << endl;
		cout << " ---O(V) \n";
		for (i = 0; i < n; i++) {
			cout << "Vertices[" << i << "] : " << Vertices[i] << endl;
		}
		cout << " ---O(E) \n";
		for (i = 0; i < m; i++) {
			cout << "AdjacencyList[" << i << "] : " << AdjacencyList[i] << endl;
		}

		cout << " Transpose Adj List " << endl;
		cout << " ---O(V) \n";
		for (i = 0; i < n; i++) {
			cout << "Vertices[" << i << "] : " << Vertices_Transpose[i] << endl;
		}
		cout << " ---O(E) \n";
		for (i = 0; i < m; i++) {
			cout << "AdjacencyList[" << i << "] : " << AdjacencyList_Transpose[i] << endl;
		}
	}

	fw_bw(n, m, Vertices, AdjacencyList, Vertices_Transpose, AdjacencyList_Transpose);

	delete(Vertices);
	delete(AdjacencyList);
	delete(Vertices_Transpose);
	delete(AdjacencyList_Transpose);
}