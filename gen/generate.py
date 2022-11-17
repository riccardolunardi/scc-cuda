import random
import concurrent.futures
import itertools
import operator

import networkx as nx
import numpy as np

from statistics import mean

CYCLES_TO_CREATE = 20
GRAPH_NAME = "./gen/sample_test"
N_ARCS_TO_REMOVE = 40
RATIO_U_NODES_TO_TOTAL_NODES = 1/5
LOWER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*2/11)
UPPER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*2/9)

def gen_cycles(_):
    cycle_length = random.randint(LOWER_BOUND_LENGTH_CYCLE, UPPER_BOUND_LENGTH_CYCLE)
    return nx.cycle_graph(cycle_length, create_using=nx.DiGraph()), cycle_length

def add_edges_to_main_graph(new_edges_with_offset):
    main_graph.add_edges_from(new_edges_with_offset)

def get_new_edges(offset, i, e):
    return ((e[0]+i)*offset, (e[1]+i)*offset)

def graphs_equal(graph1, graph2):
    return (
        graph1.adj == graph2.adj
        and graph1.nodes == graph2.nodes
        and graph1.graph == graph2.graph
    )

def add_offset_to_edges(cycle_offset):
    cycle, offset = cycle_offset
    return [(e[0]+offset, e[1]+offset) for e in cycle.edges]

if __name__ == "__main__":
    
    cycles = []
    chosen_length_cycle = []

    print("Genero i cicli...")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        cycles = executor.map(gen_cycles, range(1, CYCLES_TO_CREATE+1))
    
    cycles = list(cycles)
    cycles, offsets = zip(*cycles)
    cycles = list(zip(cycles, np.cumsum(list(offsets))))

    avg_length_cycles = int(mean([c[1] for c in cycles]))

    print(f"Cicli creati: {CYCLES_TO_CREATE}", f"Lunghezza media dei cicli creati: {avg_length_cycles}", "=================================================", sep="\n")

    print("Creazione del grafo principale...")
    main_graph = nx.fast_gnp_random_graph(int(avg_length_cycles*CYCLES_TO_CREATE/1000), p=0.001, directed=True)

    print("Unisco i cicli...")
    new_edges = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        new_edges = list(executor.map(add_offset_to_edges, cycles))
        new_edges = list(itertools.chain.from_iterable(new_edges))
        print(f"Nuovi archi calcolati {len(new_edges)}...")

        main_graph.add_edges_from(new_edges)

    print("Aggiunti archi randomici")

    print(f"Rimuovo {N_ARCS_TO_REMOVE} archi a caso...")
    main_graph.remove_edges_from(random.choices(list(main_graph.edges), k=N_ARCS_TO_REMOVE))

    print("Creazione grafo conclusa! Calcolo delle SCC...", "=================================================", sep="\n")

    sccs = [g for g in nx.strongly_connected_components(main_graph)]
    print(f"Nodi della rete: {len(main_graph.nodes)}", f"Archi della rete: {len(main_graph.edges)}", f"SCC create: {len(sccs)}", 
        f"Grandezza media delle SCC: {int(mean([len(s) for s in sccs if len(s)>2]))}",sep="\n")


    order_adj = []
    operator.itemgetter(1,2)
    for u, connected_to_u in main_graph.adjacency():
        for v in connected_to_u:
            order_adj.append((u, v))

    order_adj.sort(key=operator.itemgetter(0,1))
    
    n_nodes = len(main_graph.nodes)
    random_u_nodes = random.choices(range(n_nodes),k=int(RATIO_U_NODES_TO_TOTAL_NODES*n_nodes))

    sample_text = f"% {len(order_adj)} {len(main_graph.nodes)}\n"
    sample_text += "\n".join([f"{u} {v}" for u, v in order_adj])
    sample_text += "\n" + "\n".join(list(map(str, random_u_nodes)))
    with open(GRAPH_NAME, "w") as sample:
        sample.write(sample_text)
