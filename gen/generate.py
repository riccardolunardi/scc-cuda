from functools import partial
import random
from multiprocessing import Pool
import itertools
import operator
from itertools import repeat
import os

import networkx as nx
import numpy as np

from statistics import mean

CYCLES_TO_CREATE = 300
GRAPH_NAME = f"{os.getcwd()}/sample_test"
RANDOM_ARCS_TO_ADD = 0.99
N_ARCS_TO_REMOVE = 1500
RATIO_U_NODES_TO_TOTAL_NODES = 5/6
LOWER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*5/60)
UPPER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*59/60)

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

def keep_arc_if_node_exists(u_v, n):
    if u_v[0] > n or u_v[1] > n:
        return None
    return u_v

def remove_offset_nodes(v, offset):
    return v-offset

def remove_offset_arcs(u_v, offset):
    return (u_v[0]-offset, u_v[1]-offset)

def limit_arcs(u_v, max_node):
    if max_node >= u_v[0] and max_node >= u_v[1]:
        return u_v
    return None

if __name__ == "__main__":
    
    cycles = []
    chosen_length_cycle = []

    print("Genero i cicli...")

    cycles = map(gen_cycles, range(1, CYCLES_TO_CREATE+1))
    
    cycles = list(cycles)
    cycles, offsets = zip(*cycles)
    cycles = list(zip(cycles, np.cumsum(list(offsets))))

    node_of_every_cycle = cycles[-1][1]
    avg_length_cycles = int(node_of_every_cycle/len(cycles))

    print(f"Cicli creati: {CYCLES_TO_CREATE}", f"Lunghezza media dei cicli creati: {avg_length_cycles}", "=================================================", sep="\n")

    print("Creazione del grafo principale...")
    main_graph = nx.empty_graph(node_of_every_cycle, create_using=nx.DiGraph)

    print("Unisco i cicli...")
    new_edges = list(map(add_offset_to_edges, cycles))
    new_edges = list(itertools.chain.from_iterable(new_edges))
    
    for _ in range(int(RANDOM_ARCS_TO_ADD*node_of_every_cycle)):
        new_edges.append((random.choice(range(node_of_every_cycle)), random.choice(range(node_of_every_cycle))))

    print(f"Nuovi archi calcolati {len(new_edges)}...")

    main_graph.add_edges_from(new_edges)

    print("Aggiunti archi randomici!")
    
    print(f"Rimuovo {N_ARCS_TO_REMOVE} archi a caso...")
    main_graph.remove_edges_from(random.choices(list(main_graph.edges), k=N_ARCS_TO_REMOVE))

    print(f"Rimozione dell'offset nei ({len(main_graph.nodes)}) nodi...")
    offset = list(main_graph.nodes)[0]
    max_node_to_exists = max(list(map(remove_offset_nodes, main_graph.nodes, repeat(offset))))
    main_graph.nodes = list(range(max_node_to_exists))

    print(f"Rimozione dell'offset negli archi...")
    #edges_no_offset = list(map(remove_offset_arcs, main_graph.edges, repeat(offset)))
    with Pool() as p:
        edges_no_offset = list(p.map(partial(limit_arcs, max_node=max_node_to_exists), main_graph.edges))
    
    main_graph.edges = list(filter(lambda x: x is not None, edges_no_offset))

    final_graph = nx.DiGraph()
    final_graph.update(edges=main_graph.edges, nodes=main_graph.nodes)

    print("Creazione grafo conclusa! Calcolo delle SCC...", "=================================================", sep="\n")

    sccs = [g for g in nx.strongly_connected_components(final_graph)]
    not_only_one = list(filter(lambda x: len(x) > 1, sccs))
    print(f"Nodi della rete: {len(final_graph.nodes)}", f"Archi della rete: {len(final_graph.edges)}", f"SCC create: {len(sccs)}", 
        f"Grandezza media delle SCC: {int(mean([len(s) for s in sccs if len(s)>2]))}",sep="\n")

    order_adj = []
    for u, connected_to_u in final_graph.adjacency():
        for v in connected_to_u:
            order_adj.append((u, v))

    order_adj.sort(key=operator.itemgetter(0,1))
    
    n_nodes = len(final_graph.nodes)
    scc_in_u = itertools.chain.from_iterable(random.choices(not_only_one, k=int(0.10*len(not_only_one))))
    random_u_nodes = set(random.choices(range(n_nodes),k=int(RATIO_U_NODES_TO_TOTAL_NODES*n_nodes))).union(set(scc_in_u))

    sample_text = f"% {len(order_adj)} {len(final_graph.nodes)}\n"
    sample_text += "\n".join([f"{u} {v}" for u, v in order_adj])
    sample_text += "\n" + "\n".join(list(map(str, random_u_nodes)))
    with open(GRAPH_NAME, "w") as sample:
        sample.write(sample_text)
    