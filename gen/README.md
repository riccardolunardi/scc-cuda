# Informazioni sulla creazione di reti artificiali

## Algoritmo

L'algoritmo _molto ad alto livello_ per la creazione di grafi con un numero di componenti fortemente connesse (SCC) variabili è la seguente:

1. Generare un numero di reti cicliche di lunghezza variabile
2. Creare un grafo G che abbia un numero possibilmente variabile di archi randomici tra i nodi
3. Unire i grafi ciclici e il grafo G, in modo da ottenere un'unica grande rete H.
4. Rimuovere/Aggiungere un numero di archi/nodi variabile dal grafo H

## Implementazione

L'implementazione è avvenuta tramite il package `NetworkX` di Python. Creare grafi ciclici e grafi randomici è stato relativamente facile grazie alle funzione già predisposte.

### Note implementative

#### Sovrapposizione di archi

Unendo i vari grafi c'è il rischio di aggiungere archi e nodi già presenti, visto la numerazione per tutti comincia da 0. Questo non crea difficoltà unendo i grafi ciclici con il grafo randomico, perché vogliamo che a quest'ultimo, oltre ai suoi archi, vengano aggiunti i nodi e soprattutto gli archi dei grafi ciclici.

Il problema c'è quando si uniscono i grafi ciclici, in quanto si andrebbero a sovrapporre nodi e archi con lo stesso numero. Per evitare questo, al grafo principale aggiungiamo solo la lista di archi, che sono implementati come `tuple`, in cui il valore di ogni arco viene aggiunto un offset tale che nessun arco si sovrappone all'altro.