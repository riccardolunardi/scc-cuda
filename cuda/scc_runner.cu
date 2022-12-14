#include "../utils/file2graph.cpp"
#include "sccv4_pinned.cu"
using namespace std;

int main(unsigned int argc, char ** argv) {
    if (argc < 4) {
		cout << " Invalid Usage !! Usage is ./main.out <graph_input_file> <number_of_repetition> <profiling(0|1)>\n";
		return -1;
	}

	unsigned num_nodes, num_edges;
    unsigned * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	char * status;

    create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);

	char * og_status;
	int repeat = atoi(argv[2]);
	bool profiling = atoi(argv[3]);
	og_status = (char *) malloc(num_nodes * sizeof(char));
	memcpy(og_status, status, num_nodes);

	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		routine(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
	}

}