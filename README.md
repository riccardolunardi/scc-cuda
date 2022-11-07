# Progetto di Architetture Parallele

Richiesta del testo in data: 28 settembre 2022

Risposta con testo in data: 03 ottobre 2022

## Richiesta

Richiamo alcune **definizioni** preliminari:
  - 1 grafo diretto: G(V,E) con V insieme dei nodi, con E insieme di archi
  - 1 SCC di G: 1 qualsiasi insieme S (massimale con almeno 2 nodi) tale che per ogni coppia di nodi (x,y) in SCC esiste sempre un cammino diretto
  - Dato 1 insieme di nodi S, prec(S) = { n | n in V-S tale che esiste un arco (n,m) per almeno un qualsiasi nodo m in S }. Quindi prec(S) è l'insieme di tutti i nodi esterni all'insieme S che hanno almeno un arco connesso ad un nodo di S.

Il problema da risolvere in CUDA è il seguente: \
Dato un grafo G(V,E) e un insieme di nodi U  (un qualsiasi sottoinsieme di V),
determinare se esiste almeno una SCC  S del grafo tale che S è sottoinsieme di U
e nessun elemento di prec(S) appartiene a U.

## Vincoli da rispettare

Svolgerlo in CUDA e, opzionalmente, anche in CUDA+OpenMP.

La implementazione dovrebbe essere in grado di gestire grafi di dimensioni arbitrariamente grandi,
nei limiti della memoria disponibile sulla GPU.

Provare ad iniziare, impostando una possibile soluzione, rappresentando i dati su GPU e le principali funzioni da implementare,  
per poi riparlarne prima di iniziare a scrivere codice.

È da ricordare che la data della discussione può essere diversa dagli appelli ufficiali. \
È necessario fissare un giorno e la soluzione (codice, relazione, risultati di esempi, ...) va inviata almeno 10 giorni prima della data dell'orale.

## Domande da fare al docente

  -  ~~Il grafo è pesato? E se si, possono anche essere negativi i pesi?~~ \
  Non importa dato che per la SCC non serve sapere se è pesato o no, ma basta sapere se l'arco c'è o no
  - Quando dice "nessun elemento di prec(C) appartiene a U", in realtà vuole dire "nessun elemento di prec(S) appartiene a U"?
  
## Input
  - V = Insieme di nodi
  - E = Insieme di archi
  - U = Un sottoinsieme di V

## Output
  - Almeno 1 SCC, che chiameremo S se esiste, tale che:
    - S sia un sottoinsieme di U
    - Nessun elemento di prec(S) appartiene a U

## Come rappresentare i dati

Per aiutarci, faremo uso del [sito](https://csacademy.com/app/graph_editor/) per rappresentare i grafi visualmente partendo da una lista di adiacenze.


Si può rappresentare un grafo (in questo caso non connesso e non pesato) G(V, E) dove:

    V = {1,2,3,4,5}
    E = { (1,4), (3,5), (2,4), (4,5), (1,3) }

Con 2 rappresentazioni:
  - La matrice delle adiacenze

            [ 0 0 1 1 0 ]
            [ 0 0 0 1 0 ]
        G = [ 1 0 0 0 1 ]
            [ 1 1 0 0 1 ]
            [ 0 0 1 1 0 ]

  - Le liste delle adiacenze.
  Saranno la nostra soluzione per quanto riguarda come rappresentare i dati.

        1 = [ 3 ]   2 = [ 4 ]   3 = [ 1 ]   4 = [ 1 ]   5 = [ 3 ]
            [ 4 ]                   [ 5 ]       [ 2 ]       [ 4 ]
                                                [ 5 ]

Bisogna capire quale delle 2 rappresentazioni può essere migliore per il nostro caso in CUDA.

Useremo le liste di adiacenza.
Verranno rappresentate da 2 vettori, dove dati degli archi (u,v):

1. Il primo vettore, con lunghezza |V|, contiene i nodi u in maniera univoca.
2. Il secondo vettore, con lunghezza |E|, contiene i nodi v di tutti gli archi e sarà ordinato in base all'ordine dei nodi u nel primo vettore.

I nodi v punteranno ai rispettivi nodi u e viceversa.

Per esempio:

![Missing image data representation](img/example_representation_data.jpeg "Data representation")

## Algoritmo

Le operazioni identificate essere le più importanti sono le seguenti:
 * Generazione delle liste di adiacenza
 * Trovare le tutte SCC del grafo
 * Eliminare le SCC che non contengono nodi di U
 * Controllare che le SCC valide non abbiano nodi U che fanno parte dell'insieme prec(S)

Inizialmente si implementerà una versione che esegua queste operazioni in modo separato.

Possibili ottimizzazioni:
 * Si può evitare di evitare di calcolare le SCC in cui ci sono nodi che non fanno parte di U? _(Magari trimmandoli direttamente)_?

### Algoritmo (_apparentemente_) naive

L'algoritmo più facile è quello di effettuare le quattro operazioni in fasi sequenziali:

**Strutture dati aggiuntive**

* `removed`: vettore |V| che indica se un nodo potrebbe far parte o no di una SCC. Piuttosto di ricalcolare delle nuove liste senza un certo nodo, usiamo questo vettore di flag.
* `visited`: vettore |V| che indica se un nodo è stato visitato o no almeno da una visita forward o backward 
* `color`: Anche detto _range_, il vettore color ha lunghezza |V| e assegna ad ogni nodo un "colore" che identifica uno specifico sottografo. Due sottografi diversi hanno due colori diversi

_Opzionali_:
* `expanded`: vettore |V| che indica se un nodo è sgià stato attraversato o no

~~~

// Fase 1: Creazione delle liste di adiacenza
adj_list = generate_adj_list(G)
trasposed_adj_list = generate_trasposed_adj_list(G)

// Fase 2: Calcolo delle SCC
TRIM(G, SCC, removed)
PIVOT-GEN(G, SCC, color, removed)
do
  FWD-REACH(G, SCC, color, removed, visited, expanded)
  BWD-REACH(G, SCC, color, removed, visited, expanded)
  TRIM(G, SCC, removed)
  UPDATE(G, SCC, color, removed, visited, expanded)
  PIVOT-GEN(G, SCC, color, removed, visited)
until (nessun pivot generato)

// Fase 3: Eliminazione delle SCC che non hanno nodi che fanno parte di U
FILTER-SCC(G, SCC)

// Fase 4: Individuazione delle SCC con prec(S) ∩ U = {}
prec = PREC-GEN(G, SCC)
valid_sccs = PREC-CHECK(G, SCC, PREC)

~~~

`valid_sccs` dovrebbe essere il risultato finale