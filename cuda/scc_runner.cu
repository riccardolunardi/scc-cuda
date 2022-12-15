#include "../main.cpp"
#include "sccv1_naive.cu"
#include "sccv2_status.cu"
#include "sccv3_streams.cu"
#include "sccv4_pinned.cu"
#include "sccv5_openmp.cu"
#include "sccv6_optreach.cu"
using namespace std;

int main(unsigned int argc, char ** argv) {
    if (argc < 4) {
		cout << "Invalid Usage !! Usage is ./main.out <graph_input_file> <number_of_repetition> <profiling(0|1)>\n";
		return -1;
	}

	const int repeat = atoi(argv[2]);
	const bool profiling = atoi(argv[3]);

	// Inizializzazione di struttrure dati per la versione 1
	int num_nodes_v1, num_edges_v1;
    int * nodes_v1, * adjacency_list_v1, * nodes_transpose_v1, * adjacency_list_transpose_v1;
	bool * is_u;
	
	printf("Lettura grafo da file %s\n", argv[1]);
	printf("Lettura per versione naive...\n");
	create_graph_from_filename_v1(argv[1], num_nodes_v1, num_edges_v1, nodes_v1, adjacency_list_v1, nodes_transpose_v1, adjacency_list_transpose_v1, is_u);
	
	// Inizializzazione di struttrure dati per versione 2
	int num_nodes_v2, num_edges_v2;
    int * nodes_v2, * adjacency_list_v2, * nodes_transpose_v2, * adjacency_list_transpose_v2;
	char * og_status_v2;
	
	printf("Lettura per versione status....\n");
	create_graph_from_filename_notuns(argv[1], num_nodes_v2, num_edges_v2, nodes_v2, adjacency_list_v2, nodes_transpose_v2, adjacency_list_transpose_v2, og_status_v2);

	// Inizializzazione di struttrure dati per le versioni > 1
	printf("Lettura per versioni >= 3 e main.cpp...\n");
	unsigned num_nodes, num_edges;
    unsigned * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	char * og_status;
	
	create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, og_status);

	char * status, * status_v2;
	status = (char *) malloc(num_nodes * sizeof(char));
	status_v2 = (char *) malloc(num_nodes * sizeof(char));

	printf("Versione 0 - main.cpp\n");
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

	printf("Versione 1 - Naive\n");
	for(int i=0;i<repeat;i++){
		routine_v1(profiling, num_nodes_v1, num_edges_v1, nodes_v1, adjacency_list_v1, nodes_transpose_v1, adjacency_list_transpose_v1, is_u);
	}

	printf("Versione 2 - Status\n");
	for(int i=0;i<repeat;i++){
		memcpy(status_v2, og_status_v2, num_nodes); 
		routine_v2(profiling, num_nodes_v2, num_edges_v2, nodes_v2, adjacency_list_v2, nodes_transpose_v2, adjacency_list_transpose_v2, status_v2);
	}

	printf("Versione 3 - Streams\n");
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine_v3(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

	printf("Versione 4 - Pinned\n");
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine_v4(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

	printf("Versione 5 - OpenMP\n");
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine_v5(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

	printf("Versione 6 - Reach Ottimizzato\n");
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine_v6(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

}