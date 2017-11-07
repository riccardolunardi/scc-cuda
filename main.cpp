/*
 * main.cpp
 *
 *  Created on: 05-Nov-2017
 *      Author: divya
 */


#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#include "adj_list.cpp"
#include "reverse_adj_list.cpp"

#define degub 0

using namespace std;

__global__ void forward_closure(int *dVertices, int *dAdjList, int *subgraph, bool *visitedF, bool *terminateF, int numVertices, int numEdges)
{
	int i, ListPointer1, ListPointer2;
	
	int id = blockIdx.x*blockDim.x + threadIdx.x;
    int pivot = id+1;

	if(visitedF[id] )
	{
		ListPointer1 = dVertices[pivot-1];
		if(pivot == numVertices)
			ListPointer2 = numEdges;
		else
			ListPointer2 = dVertices[pivot];
	
		for(i = ListPointer1; i < ListPointer2; i++)
		{
			if(visitedF[i-1] == false && subgraph[pivot-1] == subgraph[i-1])
			{
				visitedF[i-1] = true;			
				terminateF = false;		
			}
		}

	}
}

void fw_bw(int n, int m, int *Vertices, int *AdjacencyList, int *Vertices_Transpose, int *AdjacencyList_Transpose)
{
	int *dVertices, *dAdjList, *dVerticesT, *dAdjListT;
	int *subgraph;
	bool *visitedF, *visitedB;
	bool *terminateF, *terminateB;

	cudaMalloc((void**)&dVertices, n*(sizeof(int)));
	cudaMalloc((void**)&dAdjList, m*(sizeof(int)));
	cudaMalloc((void**)&dVerticesT, n*(sizeof(int)));
	cudaMalloc((void**)&dAdjListT, m*(sizeof(int)));
	cudaMalloc((void**)&subgraph, n*(sizeof(int)));
	cudaMalloc((void**)&visitedF, n*(sizeof(bool)));
	cudaMalloc((void**)&visitedB, n*(sizeof(bool)));
	cudaMalloc((void**)&terminateF, 1*sizeof(bool));
	cudaMalloc((void**)&terminateB, 1*(sizeof(bool)));

	cudaMemset(subgraph, 0, n);	 
	cudaMemset(visitedF, false, n);	 
	cudaMemset(visitedB, false, n);	 
	cudaMemset(terminateF, false, 1); //pinned	 
	cudaMemset(terminateB, false, 1);	 
	
	cudaMemcpy(dVertices, Vertices, sizeof(int)*n, cudaMemcpyHosttoDevice);
	cudaMemcpy(dAdjList, AdjacencyList, sizeof(int)*m cudaMemcpyHosttoDevice);
	cudaMemcpy(dVerticesT, Vertices_Transpose, sizeof(int)*n, cudaMemcpyHosttoDevice);
	cudaMemcpy(dAdjListT, AdjacencyList_Transpose, sizeof(int)*m, cudaMemcpyHosttoDevice);

	int numBlocks, numThreadsPerBlock, pivot;
	numThreadsPerBlock = 64;
	numBlocks = n/numThreadsPerBlock + (n%numThreadsPerBlock == 0 ? 0 : 1);

	pivot = 0; 
	cudaMemset(&visitedF[pivot],true,1);
	cudaMemset(&visitedB[pivot],true,1);
	while(!terminateF)
	{
	  terminateF = true;
      forward_closure<<<numBlocks, numThreadsPerBlock>>>(pivot, dVertices, dAdjList, subgraph, visitedF, terminateF, n, m);
	}

	backward_closure<<<numBlocks, numThreadsPerBlock>>>(pivot, dVerticesT, dAdjListT, subgraph, visitedB, terminateB, n, m);

}

int main(int argc, char ** argv)
{
        if(argc != 2)
        {
                cout << " Invalid Usage !! Usage is ./a.out <graph_input_file> \n";
                return -1;
        }
        const char *filename = argv[1];

        std::string line;
        std::ifstream infile(filename);

        std::getline(infile,line);

        std::getline(infile,line);
        std::istringstream iss(line);

        char c;
        int m, n, x, i;

        iss >> c;
        iss >> m;
        iss >> n;
        iss >> x;

        std::cout << "number of vertices " << n << endl;
        std::cout << "number of edges " << m << endl;

        int *Vertices;
        int *AdjacencyList;

        int *Vertices_Transpose;
        int *AdjacencyList_Transpose;

        Vertices = new int[n];
        AdjacencyList = new int[m];

        Vertices_Transpose = new int[n];
        AdjacencyList_Transpose = new int[m];

        infile.close();

        adj_list(filename, Vertices, AdjacencyList);
        reverse_adj_list(filename, Vertices_Transpose, AdjacencyList_Transpose);

	if(debug)
	{
       	cout << " Adj List " << endl;
        cout << " ---O(V) \n";
        for(i = 0; i < n; i++)
        {
	        cout << "Vertices["<< i << "] : " << Vertices[i] << endl;
        }
        cout << " ---O(E) \n";
        for(i = 0; i < m; i++)
        {
        	cout << "AdjacencyList[" << i << "] : " << AdjacencyList[i] << endl;
        }

        cout << " Transpose Adj List " << endl;
        cout << " ---O(V) \n";
        for(i = 0; i < n; i++)
        {
    	    cout << "Vertices["<< i << "] : " << Vertices_Transpose[i] << endl;
        }
        cout << " ---O(E) \n";
        for(i = 0; i < m; i++)
        {
        	cout << "AdjacencyList[" << i << "] : " << AdjacencyList_Transpose[i] << endl;
        }
	}

	fw_bw(n, m, Vertices, AdjacencyList, Vertices_Transpose, AdjacencyList_Transpose);	

	delete(Vertices);
	delete(AdjacencyList);
	delete(Vertices_Transpose);
	delete(AdjacencyList_Transpose);
}
