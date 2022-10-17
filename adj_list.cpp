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

/*
Informazioni base da ricordare:
 - n = numero di vertici del grafo
 - m = numero di archi del grafo
 - iss è una funzione che apparentemente, data una stringa, la splitta in base agli spazi (" ") e con
   l'operatore >> restituisce i valori uno alla volta, in ordine.
*/

void adj_list(const char *filename, int *Vertices, int *AdjacencyList)
{
	// La prima riga viene scartata come nel main
	std::string line;
	std::ifstream infile(filename);
	std::getline(infile, line);

	// Nella seconda vengono prese le principali informazioni
	std::getline(infile, line);
	std::istringstream iss(line);

	// Questa è inutile
	char c;
	iss >> c;

	// Inizializzazione di tutte le variabili intere: archi, nodi, variabile di debuff, u e v
	// - u rappresenterà il numero del vertice da cui parte l'arco che si sta elaborando
	// - v rapprenseterà il numbero del vertice a cui l'arco che si sta elaborando punta
	int m, n, x, u, v;
	int i;

	iss >> m;
	iss >> n;
	iss >> x;

	DEBUG_MSG("---- Obtaining the adjacency list ----");
	DEBUG_MSG("Number of edges: " << m);
	DEBUG_MSG("Number of vertices: " << n);

	// uedges è il numero di arco che si sta elaborando nel ciclo corrente
	// int uedges = 0;
	int oldu = 0;

	//-------------------------------------------------------------------------
	//-------------------------Filling O(V) list-------------------------------
	//-------------------------------------------------------------------------
	/*
	Questa parte di codice è un po' confusa e suppone che gli archi nel file ("u v")
	siano scritti in modo tale che "u" vengano riportati dal valore di "u" più piccolo al
	valore di "u" più grande, in ordine crescente.
	Il risultato è lista una lista di lunghezza |V|, a cui ogni nodo (identificato da l'indice) viene
	associato il numero di nodi uscenti (out degree)
	*/
	while (std::getline(infile, line))
	{

		std::istringstream iss(line);

		iss >> u;
		iss >> v;

		Vertices[u - 1] += 1;

		//Debuffing, si legge finché non c'è niente
		while (iss >> x)
		{
		}
	}

	if(debug){
		cout << "Counting how many arcs a node has (temporary)" << endl;
		for (i = 0; i < n; i++)
		{
			cout << "Vertices[" << i << "]: " << Vertices[i] << endl;
		}
		cout << endl;
	}

	/*
	In questa sezione, quello che prima era l'out degree di un singolo nodo, ora diventa l'out degree
	del singolo nodo + la somma di tutti gli out degree dei nodi con indice minore del suo. Questo è fatto
	in modo tale da permettere alla struttura dati di poter accedere (tramite questo valore appena calcolato)
	al giusto indice della seconda lista.
	*/
	int temp1, temp2 = 0;
	for (i = 1; i < n; i++)
	{
		temp1 = Vertices[i];
		Vertices[i] = Vertices[i - 1] + temp2;
		temp2 = temp1;
	}

	//Il primo indice sarà sempre uguale a 0
	Vertices[0] = 0;

	if(debug){
		cout << "----O(V) list----" << endl;
		for (i = 0; i < n; i++){
			cout << "Vertices[" << i << "] : " << Vertices[i] << endl;
		}
		cout << endl << endl;
	}
	
	infile.close();

	//-------------------------------------------------------------------------
	//-------------------------Filling O(E) list-------------------------------
	//-------------------------------------------------------------------------

	infile.open(filename);
	// Si scartano le prime due righe
	std::getline(infile, line);
	std::getline(infile, line);

	int idx = 0;
	oldu = -1;

	while (std::getline(infile, line)) // Finché nel file ci sono righe da leggere...
	{
		std::istringstream iss(line);

		// Lettura degli archi
		iss >> u;
		iss >> v;
		// È vero che l'ultimo nodo elaborato dal ciclo è uguale a quello corrente?
		if (oldu == u - 1)
		{
			// Se ci capita qui, vuol dire che si sta lavorando sempre con lo stesso nodo "u", quindi
			// di aggiungono in successione i nodi puntati da "u"
			AdjacencyList[idx++] = v;

			// std::cout << "Entering "<< u <<" to "<< v << " at "<< idx-1;
		}
		else
		{
			// idx acquisisce la posizione dalla quale verranno salavati i nodi puntati da u
			idx = Vertices[u - 1];
			AdjacencyList[idx++] = v;
			oldu = u - 1;
			// std::cout << "Entering " << u <<" to "<< v << " at "<< idx-1;
		}

		// Debuffing come prima
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
		cout << "-----------------------------------------" << endl;
	}

	infile.close();
}