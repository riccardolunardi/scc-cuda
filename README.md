# Progetto di Architetture Parallele

Ricevuto il testo in data: 03 ottobre 2022

## Istruzioni per l'uso

Per la compilazione e l'esecuzione del programma, si utilizzano i seguenti comandi presenti nel MAKE file:
- `make cpp-compile` per compilare la versione serializzata dell'algoritmo
- `make cpp-run` per eseguire la versione serializzata dell'algoritmo
- `make cuda-compile` per compilare il runner di CUDA
- `make cuda-run` per eseguire il runner di CUDA

Per utilizzare i grafi discussi nella relazione, utilizzare il link fornito per scaricarli. Una volta scompattati, inserire il percorso all'interno del MakeFile nel comando corretto.

### Attenzione

Il runner (`scc_runner.cu`) esegue tutte le versioni dell'algoritmo, utilizzando il grafo fornito in input. Per eseguire solo le versioni desiderate, si commenti la porzione di codice apposita e si ricompili.

### Esempi di comandi

#### Compilazione della versione serializzata

```g++ -std=c++11 standalone.cpp -o ./build/standalone.exe```

#### Compilazione della versione parallela

```nvcc -Xcompiler /openmp -DDEBUG_FINAL=1 -DOMP_MIN_NODES=100000 -DWARMUP=0 ./cuda/scc_runner.cu -o ./build/scc.exe```

- **DEBUG_FINAL**: 0 per non far printare niente a video, 1 per visualizzare il risultato
- **OMP_MIN_NODES**: Numero di nodi necessario per attivare l'uso di OpenMP
- **WARMUP**: parametro usato per la raccolta dei tempi di esecuzione. Corrisponde al numero di esecuzioni successive di cui non viene salvata la durata.

#### Esecuzione della versione serializzata

```./build/standalone.exe F:/network-benchmark/final/twitter.txt```

Il parametro è la rete da fornire in input

#### Esecuzione della versione parallela

```./build/scc.exe F:/network-benchmark/final/twitter.txt 1```

I parametri sono rispettivamente:
- Rete da fornire in input
- Numero di esecuzioni da fare per ogni versione (senza considerare il parametro WARMUP)

## Richiesta del problema

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
