#include "../main.cpp"
#include "sccv1_naive.cu"
#include "sccv2_status.cu"
#include "sccv3_streams.cu"
#include "sccv4_pinned.cu"
#include "sccv5_openmp.cu"
#include "sccv6_optreach.cu"
#include "sccv7_optreach.cu"
#include <chrono>
#include <vector>
using namespace std;

#define WARMUP 3

double calculateStandardDeviation(double mean, int n, double numbers[]) {
	// Calculate the sum of squared differences
	double sumSquaredDifferences = 0;
	for (int i = 0; i < n; i++) {
		double difference = numbers[i] - mean;
		sumSquaredDifferences += difference * difference;
	}

	// Calculate the standard deviation
	return sqrt(sumSquaredDifferences / n);
}

vector<double> common_routine(void (*routine_runner)(const bool, unsigned int, unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*, char * ), const bool profiling, const unsigned int num_nodes, const unsigned int num_edges, unsigned int * nodes, unsigned int * adjacency_list, unsigned int * nodes_transpose, unsigned int * adjacency_list_transpose, char * status, char * og_status, const int repeat) {
	// Call the function passed as an argument
	vector<double> executionTimes;

	for(short i=0; i<repeat + WARMUP; i++){
		memcpy(status, og_status, num_nodes);

		auto start = chrono::high_resolution_clock::now();
		routine_runner(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
		auto end = chrono::high_resolution_clock::now();

		if (WARMUP < i) {
			executionTimes.push_back(chrono::duration<double, milli>(end - start).count());
		}
	}

	return executionTimes;
}

void print_benchmark(const vector<double> executionTimes) {
	// Calculate the mean
	double sum = 0;
	for (double t : executionTimes) {
		sum += t;
	}
	double mean = sum / executionTimes.size();

	// Calculate the sum of squared differences
	double sumSquaredDifferences = 0;
	for (double t : executionTimes) {
		double difference = t - mean;
		sumSquaredDifferences += difference * difference;
	}

	// Calculate the standard deviation
	double standardDeviation = sqrt(sumSquaredDifferences / executionTimes.size());

  	cout << "Total elapsed time: " << sum << "ms" << endl;
	cout << "Average elapsed time: " << mean << "ms" << endl;
	cout << "Standard deviation: " << standardDeviation << "ms" << endl;
}

int main(unsigned int argc, char ** argv) {
    if (argc < 4) {
		cout << "Invalid Usage !! Usage is ./main.out <graph_input_file> <number_of_repetition> <profiling(0|1)>\n";
		return -1;
	}

	const int repeat = atoi(argv[2]);
	const bool profiling = atoi(argv[3]);

	// Inizializzazione di struttrure dati per le versioni > 1
	unsigned num_nodes, num_edges;
    unsigned * nodes, * adjacency_list, * nodes_transpose, * adjacency_list_transpose;
	char * og_status;
	
	printf("Lettura del file %s...\n", argv[1]);
	create_graph_from_filename(argv[1], num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, og_status);

	printf("Number of nodes: %d...\n", num_nodes);
	
	// Inizializzazione di struttrure dati per la versione 1
	int num_nodes_v1, num_edges_v1;
    int * nodes_v1, * adjacency_list_v1, * nodes_transpose_v1, * adjacency_list_transpose_v1;
	bool * is_u;
	
	num_nodes_v1 = (int)num_nodes;
	num_edges_v1 = (int)num_edges;
	nodes_v1 = (int *) malloc((num_nodes_v1+1) * sizeof(int));
	adjacency_list_v1 = (int *) malloc(num_edges_v1 * sizeof(int));
	nodes_transpose_v1 = (int *) malloc((num_nodes_v1+1) * sizeof(int));
	adjacency_list_transpose_v1 = (int *) malloc(num_edges_v1 * sizeof(int));
	is_u = (bool *) malloc(num_nodes_v1 * sizeof(bool));

	for(int i=0; i<num_nodes; i++){
		nodes_v1[i] = (int)nodes[i];
		nodes_transpose_v1[i] = (int)nodes_transpose[i];
		is_u[i] = get_is_u(og_status[i]);
	}

	for(int i=0; i<num_edges; i++){
		adjacency_list_v1[i] = (int)adjacency_list[i];
		adjacency_list_transpose_v1[i] = (int)adjacency_list_transpose[i];
	}

	char * status;
	status = (char *) malloc(num_nodes * sizeof(char));

 	printf("Versione 0 - main.cpp\n");
	vector<double> executionTimes;
	for(int i=0;i<repeat;i++){
		memcpy(status, og_status, num_nodes);
		
		auto start = chrono::high_resolution_clock::now();
		routine(profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status);
		auto end = chrono::high_resolution_clock::now();
		if (WARMUP < i) {
			executionTimes.push_back(chrono::duration<double, milli>(end - start).count());
		}
	}

  	print_benchmark(executionTimes);

	printf("Versione 1 - Naive\n");
	executionTimes.clear();
	for(int i=0;i<repeat;i++){
		auto start = chrono::high_resolution_clock::now();
		routine_v1(profiling, num_nodes_v1, num_edges_v1, nodes_v1, adjacency_list_v1, nodes_transpose_v1, adjacency_list_transpose_v1, is_u);
		auto end = chrono::high_resolution_clock::now();

		if (WARMUP < i) {
			executionTimes.push_back(chrono::duration<double, milli>(end - start).count());
		}
	}

	print_benchmark(executionTimes);

	printf("Versione 2 - Status\n");
	print_benchmark(common_routine(routine_v2, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));
	
	printf("Versione 3 - Streams\n");
	print_benchmark(common_routine(routine_v3, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));
	
	printf("Versione 4 - Pinned\n");
	print_benchmark(common_routine(routine_v4, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));
	
	printf("Versione 5 - OpenMP\n");
	print_benchmark(common_routine(routine_v5, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));
	
	printf("Versione 6 - Reach Ottimizzato\n");
	print_benchmark(common_routine(routine_v6, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));

	printf("Versione 7 - Reach Ottimizzato + status unico\n");
	print_benchmark(common_routine(routine_v7, profiling, num_nodes, num_edges, nodes, adjacency_list, nodes_transpose, adjacency_list_transpose, status, og_status, repeat));
}