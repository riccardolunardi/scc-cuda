#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#define debug true
using namespace std;

#define DEBUG_MSG(str){                     	\
    if(debug)                            		\
        std::cout << str<< std::endl;         	\
}

void reverse_adj_list(const char *filename, int *nodes, int *adjacency_list) {
	std::string line;
	std::ifstream infile(filename);

	int *edges_count_per_node;

	std::getline(infile, line);
	std::getline(infile, line);
	std::istringstream iss(line);

	char percentage_sign;
	iss >> percentage_sign;

	int num_edges, num_nodes, edge_weight, first_node_of_edge, second_node_of_edge;
	int i;

	iss >> num_edges;
	iss >> num_nodes;
	iss >> edge_weight;

	DEBUG_MSG("---- Obtaining the reversed adjacency list ----");
	DEBUG_MSG("Number of edges: " << num_edges);
	DEBUG_MSG("Number of nodes: " << num_nodes);

	edges_count_per_node = new int[num_nodes];

	//-------------------------------------------------------------------------
	//-------------------------Filling O(V) list-------------------------------
	//-------------------------------------------------------------------------
	while (std::getline(infile, line)) {
		std::istringstream iss(line);

		iss >> first_node_of_edge;
		iss >> second_node_of_edge;

		nodes[second_node_of_edge - 1] += 1;

		while (iss >> edge_weight) {}
	}

	int temp1, temp2 = 0;
	for (i = 1; i < num_nodes; i++) {
		temp1 = nodes[i];
		nodes[i] = nodes[i - 1] + temp2;
		temp2 = temp1;
	}

	DEBUG_MSG("----O(V) list----");

	for (i = 0; i < num_nodes; i++) {
		if (i == 0) {
			nodes[i] = 0;
        }

		DEBUG_MSG("nodes[" << i << "] : " << nodes[i]);
		edges_count_per_node[i] = nodes[i];
	}

	DEBUG_MSG("\n");
	infile.close();

	//-------------------------------------------------------------------------
	//-------------------------Filling O(E) list-------------------------------
	//-------------------------------------------------------------------------

	infile.open(filename);
	std::getline(infile, line);
	std::getline(infile, line);

	int idx = 0;
	int old_second_node_of_edge = -1;

	// qui non ci ero mai arrivato
	while (std::getline(infile, line)) {
		std::istringstream iss(line);

		iss >> first_node_of_edge;
		iss >> second_node_of_edge;

		old_second_node_of_edge = second_node_of_edge;
		idx = edges_count_per_node[old_second_node_of_edge];

		adjacency_list[idx] = first_node_of_edge;
		edges_count_per_node[old_second_node_of_edge] += 1;

		while (iss >> edge_weight) {}
	}

	if (debug) {
		cout << "----O(E) list----\n";
		for (i = 0; i < num_edges; i++) {
			cout << "adjacency_list[" << i << "] : " << adjacency_list[i] << endl;
		}
		cout << endl << endl;
	}
	infile.close();
	delete (edges_count_per_node);
}