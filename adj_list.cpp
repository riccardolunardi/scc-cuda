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
 - n = numero di nodi del grafo
 - m = numero di archi del grafo
 - iss è una funzione che apparentemente, data una stringa, la splitta in base agli spazi (" ") e con
   l'operatore >> restituisce i valori uno alla volta, in ordine.
*/

void adj_list(const char *filename, int *nodes, int *adjacency_list) {
	// La prima riga viene scartata come nel main
	std::string line;
	std::ifstream infile(filename);
	std::getline(infile, line);

	// Nella seconda vengono prese le principali informazioni
	std::getline(infile, line);
	std::istringstream iss(line);

	// Questa è inutile
	char percentage_sign;
	iss >> percentage_sign;

	// Inizializzazione di tutte le variabili intere: archi, nodi, variabile di debuff, u e v
	// - first_node_of_edge rappresenterà il numero del nodo da cui parte l'arco che si sta elaborando
	// - v rapprenseterà il numbero del nodo a cui l'arco che si sta elaborando punta
	int num_edges, num_nodes, weight, first_node_of_edge, second_node_of_edge;
	int i;

	iss >> num_edges;
	iss >> num_nodes;
	iss >> weight;

	DEBUG_MSG("---- Obtaining the adjacency list ----");
	DEBUG_MSG("Number of edges: " << num_edges);
	DEBUG_MSG("Number of nodes: " << num_nodes);

	//-------------------------------------------------------------------------
	//-------------------------Filling O(V) list-------------------------------
	//-------------------------------------------------------------------------
	/*
	Questa parte di codice è un po' confusa e suppone che gli archi nel file ("first_node_of_edge second_node_of_edge")
	siano scritti in modo tale che "first_node_of_edge" vengano riportati dal valore di "first_node_of_edge" più piccolo al
	valore di "first_node_of_edge" più grande, in ordine crescente.
	Il risultato è lista una lista di lunghezza |V|, a cui ogni nodo (identificato da l'indice) viene
	associato il numero di nodi uscenti (out degree)
	*/
	while (std::getline(infile, line)) {
		std::istringstream iss(line);

		iss >> first_node_of_edge;
		iss >> second_node_of_edge;

		nodes[first_node_of_edge - 1] += 1;

		//Debuffing, si legge finché non c'è niente
		while (iss >> weight) {}
	}

	if(debug) {
		cout << "Counting how many arcs a node has (temporary)" << endl;
		for (i = 0; i < num_nodes; i++) {
			cout << "nodes[" << i << "]: " << nodes[i] << endl;
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
	for (i = 1; i < num_nodes; i++) {
		temp1 = nodes[i];
		nodes[i] = nodes[i - 1] + temp2;
		temp2 = temp1;
	}

	//Il primo indice sarà sempre uguale a 0
	nodes[0] = 0;

	if(debug) {
		cout << "----O(V) list----" << endl;
		for (i = 0; i < num_nodes; i++){
			cout << "nodes[" << i << "] : " << nodes[i] << endl;
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
	int old_first_node_of_edge = -1;

    // Finché nel file ci sono righe da leggere... 
	while (std::getline(infile, line)) {
		std::istringstream iss(line);

		// Lettura degli archi
		iss >> first_node_of_edge;
		iss >> second_node_of_edge;
		// È vero che l'ultimo nodo elaborato dal ciclo è uguale a quello corrente?
		if (old_first_node_of_edge == first_node_of_edge - 1) {
			// Se ci capita qui, vuol dire che si sta lavorando sempre con lo stesso nodo "first_node_of_edge", quindi
			// di aggiungono in successione i nodi puntati da "first_node_of_edge"
			adjacency_list[idx++] = second_node_of_edge;

			// std::cout << "Entering "<< first_node_of_edge <<" to "<< second_node_of_edge << " at "<< idx-1;
		} else {
			// idx acquisisce la posizione dalla quale verranno salavati i nodi puntati da first_node_of_edge
			idx = nodes[first_node_of_edge - 1];
			adjacency_list[idx++] = second_node_of_edge;
			old_first_node_of_edge = first_node_of_edge - 1;
			// std::cout << "Entering " << first_node_of_edge <<" to "<< second_node_of_edge << " at "<< idx-1;
		}

		// Debuffing come prima
		while (iss >> weight) {}
	}

	if (debug) {
		cout << "----O(E) list----\n";
		for (i = 0; i < num_edges; i++) {
			cout << "adjacency_list[" << i << "] : " << adjacency_list[i] << endl;
		}
		cout << "-----------------------------------------" << endl;
	}

	infile.close();
}