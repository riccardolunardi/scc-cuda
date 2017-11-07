/*
 * reverse_adj_list.cpp
 *
 *  Created on: 05-Nov-2017
 *      Author: divya
 */

#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#define debug 0

using namespace std;

// n = num of vertices, m = num of edges

//int main(int argc, char **argv)
void reverse_adj_list(const char *filename, int *Vertices, int *AdjacencyList)

{
//	const char *filename=argv[1];
	std::string line;
	std::ifstream infile(filename);

//	int *Vertices;
//	int *AdjacencyList;
	int *edges_count_per_vertex;

	std::getline(infile,line);

	std::getline(infile,line);
	std::istringstream iss(line);

	char c;
	iss>>c;

	int m,n,x,u,v;
	int i;

	iss>>m;
	iss>>n;
	iss>>x;

	if(debug)
	{
		std::cout << "number of edges=" <<m<< '\n';
		std::cout << "number of vertices = " <<n<< '\n';
	}
//	Vertices = new int[n];
//	AdjacencyList = new int[m];
	edges_count_per_vertex = new int[n];

	int oldu = -1;
//-------------------------------------------------------------------------
//-------------------------Filling O(V) list-------------------------------
//-------------------------------------------------------------------------
	while (std::getline(infile, line))  // this does the checking!
	{

		std::istringstream iss(line);

		iss>>v;
		iss>>u;

		Vertices[u-1] += 1;

		while (iss >> x)
		{
		}
	}

	int temp1, temp2=0;
	for(i=1;i<n;i++)
	{
		temp1 = Vertices[i];
		Vertices[i] = Vertices[i-1]+temp2;
		temp2 = temp1;
//		cout<<i<<" : "<<Vertices[i]<<"\n";
	}

	if(debug)
		cout<<"----O(V) list----\n";

	for(i=0;i<n;i++)
	{
		if(i == 0)
			Vertices[i] = 0;
//		else if(Vertices[i] == Vertices[i-1])
//			Vertices[i-1] = -1;
		if(debug)
			cout<<"Vertices["<<i<<"] : "<<Vertices[i]<<"\n";
		edges_count_per_vertex[i]=Vertices[i];
	}
	if(debug)
		cout<<"----O(V) list----\n";
	infile.close();

//-------------------------------------------------------------------------
//-------------------------Filling O(E) list-------------------------------
//-------------------------------------------------------------------------

	infile.open(filename);
	std::getline(infile,line);
	std::getline(infile,line);

	int idx = 0;
	oldu = -1;

	while (std::getline(infile, line))  // this does the checking!
	{
		std::istringstream iss(line);

		iss>>v;
		iss>>u;

		oldu = u-1;
		idx = edges_count_per_vertex[oldu];

//		cout << " count["<<oldu<<"] : "<< edges_count_per_vertex[oldu] << "idx : " << idx << endl;
		AdjacencyList[idx] = v;
		edges_count_per_vertex[oldu] += 1;

//		std::cout << "else AdjacencyList["<<idx<<"] : "<<v << '\n';

		while (iss >> x)
		{
		}
	}

	if(debug)
	{
		cout<<"----O(E) list----\n";
		for(i=0;i<m;i++)
		{
			cout<<"AdjacencyList["<<i<<"] : " << AdjacencyList[i]<<endl;
		}
		cout<<"----O(E) list----\n";
	}
	infile.close();
	delete(edges_count_per_vertex);
}
