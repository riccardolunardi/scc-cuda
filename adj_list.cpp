/*
 * adj_list.cpp
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
void adj_list(const char *filename, int *Vertices, int *AdjacencyList)
{
//	const char *filename=argv[1];
	std::string line;
	std::ifstream infile(filename);

//	int *Vertices;
//	int *AdjacencyList;

	std::getline(infile,line);

	std::getline(infile,line);
	std::istringstream iss(line);

	char c;
	iss>>c;

	int m,n,x, u, v;
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

	int uedges = -1;

	int oldu = -1;
//-------------------------------------------------------------------------
//-------------------------Filling O(V) list-------------------------------
//-------------------------------------------------------------------------
	while (std::getline(infile, line))  // this does the checking!
	{

		std::istringstream iss(line);

		iss>>u;
//		cout<< "oldu : "<< oldu<<", u-1 : "<<u-1<<"\n";
		if(oldu == (u-1))
		{
			uedges++;
		}
		else
		{
			if(oldu >= 0)
				Vertices[oldu] = uedges;
//			cout<<" Vertices["<<oldu<<"] : " << Vertices[oldu]<<"\n";
			uedges = 1;
			oldu = u-1;
		}

		while (iss >> x)
		{
		}
	}

	Vertices[oldu]=uedges;
//	cout<<" Vertices["<<oldu<<"] : " << Vertices[oldu]<<"\n";

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

		iss>>u;
		iss>>v;

		if(oldu == u-1)
		{
			AdjacencyList[idx++] = v;
			// std::cout << "entering "<<u <<" to "<<v << " at "<<idx-1 << '\n';
		}
		else
		{
			idx = Vertices[u-1];
			AdjacencyList[idx++] = v;
			oldu = u-1;
			// std::cout << "entering "<<u <<" to "<<v << " at "<<idx-1 << '\n';
		}

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
}

