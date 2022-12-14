#include "sccv3_streams.cu"
using namespace std;

int main(unsigned int argc, char ** argv) {
    if (argc < 4) {
		cout << "Invalid Usage !! Usage is ./main.out <graph_input_file> <number_of_repetition> <profiling(0|1)>\n";
		return -1;
	}

	const int repeat = atoi(argv[2]);
	const bool profiling = atoi(argv[3]);

	// Inizializzazione di struttrure dati per la versione 1
	/* int num_nodes_v1, num_edges_v1;
    int * nodes_v1, * adjacency_list_v1, * nodes_transpose_v1, * adjacency_list_transpose_v1;
	bool * is_u;  */
	
	//create_graph_from_filename_v1(argv[1], num_nodes_v1, num_edges_v1, nodes_v1, adjacency_list_v1, nodes_transpose_v1, adjacency_list_transpose_v1, is_u);
	
	// Inizializzazione di struttrure dati per versione 2
	/* int num_nodes_v2, num_edges_v2;
    int * nodes_v2, * adjacency_list_v2, * nodes_transpose_v2, * adjacency_list_transpose_v2; */
	char * status;
	
	//create_graph_from_filename_notuns(argv[1], num_nodes_v2, num_edges_v2, nodes_v2, adjacency_list_v2, nodes_transpose_v2, adjacency_list_transpose_v2, status);

	// Inizializzazione di struttrure dati per le versioni > 1
	unsigned num_nodes, num_edges;
    unsigned * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	//char * status;
	
	create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);

	char * og_status;
	og_status = (char *) malloc(num_nodes * sizeof(char));
	memcpy(og_status, status, num_nodes); 

	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		//routine_v1(profiling, num_nodes_v1, num_edges_v1, nodes_v1, adjacency_list_v1, nodes_transpose_v1, adjacency_list_transpose_v1, is_u);
		//routine_v2(profiling, num_nodes_v2, num_edges_v2, nodes_v2, adjacency_list_v2, nodes_transpose_v2, adjacency_list_transpose_v2, status);
		routine(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

}