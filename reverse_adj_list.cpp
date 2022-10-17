#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#define debug 1
using namespace std;

#ifdef debug
#define DEBUG_MSG(str)                 \
	do                                 \
	{                                  \
		std::cout << str << std::endl; \
	} while (false)
#else
#define DEBUG_MSG(str) \
	do                 \
	{                  \
	} while (false)
#endif

void reverse_adj_list(const char *filename, int *Vertices, int *AdjacencyList)
{
	std::string line;
	std::ifstream infile(filename);

	int *edges_count_per_vertex;

	std::getline(infile, line);
	std::getline(infile, line);
	std::istringstream iss(line);

	char c;
	iss >> c;

	int m, n, x, u, v;
	int i;

	iss >> m;
	iss >> n;
	iss >> x;

	DEBUG_MSG("---- Obtaining the reversed adjacency list ----");
	DEBUG_MSG("Number of edges: " << m);
	DEBUG_MSG("Number of vertices: " << n);

	edges_count_per_vertex = new int[n];

	int oldu = -1;
	//-------------------------------------------------------------------------
	//-------------------------Filling O(V) list-------------------------------
	//-------------------------------------------------------------------------
	while (std::getline(infile, line))
	{

		std::istringstream iss(line);

		iss >> v;
		iss >> u;

		Vertices[u - 1] += 1;

		while (iss >> x)
		{
		}
	}

	int temp1, temp2 = 0;
	for (i = 1; i < n; i++)
	{
		temp1 = Vertices[i];
		Vertices[i] = Vertices[i - 1] + temp2;
		temp2 = temp1;
	}

	DEBUG_MSG("----O(V) list----");

	for (i = 0; i < n; i++)
	{
		if (i == 0)
			Vertices[i] = 0;

		DEBUG_MSG("Vertices[" << i << "] : " << Vertices[i]);
		edges_count_per_vertex[i] = Vertices[i];
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
	oldu = -1;

	// qui non ci ero mai arrivato
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);

		iss >> v;
		iss >> u;

		oldu = u - 1;
		idx = edges_count_per_vertex[oldu];

		AdjacencyList[idx] = v;
		edges_count_per_vertex[oldu] += 1;

		while (iss >> x)
		{
		}
	}

	if (debug)
	{
		cout << "----O(E) list----\n";
		for (i = 0; i < m; i++)
		{
			cout << "AdjacencyList[" << i << "] : " << AdjacencyList[i] << endl;
		}
		cout << endl << endl;
	}
	infile.close();
	delete (edges_count_per_vertex);
}