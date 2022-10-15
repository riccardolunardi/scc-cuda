#include <sstream>
#include <fstream>
#include <string>
#include <iostream>

#define debug 0
using namespace std;

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
	//  u rappresenterà il numero del vertice da cui parte l'arco che si sta elaborando
	//  v rapprenseterà il numbero del vertice a cui l'arco che si sta elaborando punta
	int m, n, x, u, v;
	int i;

	iss >> m;
	iss >> n;
	iss >> x;

	if (debug)
	{
		std::cout << "Number of edges: " << m << '\n';
		std::cout << "Number of vertices: " << n << '\n';
	}

	// uedges è il numero di arco che si sta elaborando nel ciclo corrente
	int uedges = -1;
	int oldu = -1;

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
	while (std::getline(infile, line)) // Finché nel file ci sono righe da leggere
	{
		// Lettura della riga
		std::istringstream iss(line);
		iss >> u;
		// cout<< "oldu : "<< oldu<<", u-1 : "<<u-1<<"\n";

		// Se il nodo "u" che si sta leggendo è lo stesso della riga prima...
		if (oldu == (u - 1))
		{
			// Allora somma 1 ai suoi archi uscenti
			uedges++;
		}
		else
		{
			// Se questo non è la prima iterazione (oldu è minore di 0 solo alla prima iterazione del ciclo)
			if (oldu >= 0)
			{
				// Associa al nodo il numero di archi
				Vertices[oldu] = uedges;
				// cout<<" Vertices["<<oldu<<"] : " << Vertices[oldu]<<"\n";
			}

			// Reinizializza le variabili, ovvero assegnando il numero di archi base (1) e
			// impostando come oldu il nodo appena elaborato
			uedges = 1;
			oldu = u - 1;
		}

		// x è una variabile di debuff, ovvero scarta tutti gli input dopo "u v", in modo da finire la lettura della riga
		while (iss >> x)
		{
		}
	}

	// Ultimo assegnamento
	Vertices[oldu] = uedges;
	// cout<<" Vertices["<<oldu<<"] : " << Vertices[oldu]<< "\n";

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
		// cout<<i<<" : "<<Vertices[i]<<"\n";
	}

	// TODO: Check the next section of code can be replaced with this:
	/*
	Vertices[0] = 0;
	if (debug){
		cout << "----O(V) list----\n";
		for (i = 0; i < n; i++){
			cout << "Vertices[" << i << "] : " << Vertices[i] << "\n";
		}
		cout << "----O(V) list----\n";
	}
	*/

	if (debug)
		cout << "----O(V) list----\n";

	for (i = 0; i < n; i++)
	{
		if (i == 0)
			Vertices[i] = 0;
		// else if(Vertices[i] == Vertices[i-1])
		//	Vertices[i-1] = -1;

		if (debug)
			cout << "Vertices[" << i << "] : " << Vertices[i] << "\n";
	}

	if (debug)
		cout << "----O(V) list----\n";

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

			// std::cout << "Entering "<< u <<" to "<< v << " at "<< idx-1 << '\n';
		}
		else
		{
			// idx acquisisce la posizione dalla quale verranno salavati i nodi puntati da u
			idx = Vertices[u - 1];
			AdjacencyList[idx++] = v;
			oldu = u - 1;
			// std::cout << "Entering " << u <<" to "<< v << " at "<< idx-1 << '\n';
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
		cout << "----O(E) list----\n";
	}

	infile.close();
}