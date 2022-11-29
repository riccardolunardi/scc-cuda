import networkx as nx
import re
import itertools

def get_prec(G, scc):
    """
    Dato una scc e il grafo G di appartenenza, get_prec calcola l'insieme prec(scc). L'insieme prec include tutti i nodi che puntano alla scc
    :param G: grafo contienente scc
    :param scc: la scc di cui si vuole calcolare l'insieme prec
    :returns: prec(scc)
    """
    prec_with_scc = [list([i[0] for i in G.in_edges(n)]) for n in scc]
    prec_with_scc = set(itertools.chain.from_iterable(prec_with_scc))    
    return prec_with_scc - scc

nodes = set()
edges = []
u_nodes = []

# Lettura del file contenete il grafo
with open("../samples/mid_tests/sample_test_scc_fewu", "r") as main_graph_file:
    for line in main_graph_file.readlines():
        if "%" in line:
            pass
        elif re.findall(r"\w+ \w+", line):
            u, v = line.split(" ")
            u, v = int(u), int(v)
            nodes.add(u)
            edges.append((u, v))
        else:
            u_nodes.append(int(line.strip()))

# Creazione del grafo dati nodi e archi appena ottenuti dal file di test
main_graph = nx.DiGraph()
main_graph.add_nodes_from(nodes)
main_graph.add_edges_from(edges)

# Collezione e filtraggio delle SCCs
# Verrano considerate solo:
# - le SCCs non triviali verranno considerate (ovvero quelle con più di un elemento)
# - le SCCs che contengono solo nodi di U: altrimenti la SCC non è valida, secondo la richiesta del progetto
real_sccs = []
every_sccs = list(nx.strongly_connected_components(main_graph))
for scc in every_sccs:
    if len(scc) > 1:
        if all([n in u_nodes for n in scc]):
            real_sccs.append(scc)

# Qui vengono scartate le SCC che contengono nodi U nell'insieme prec(SCC)
final_sccs = []
for scc in real_sccs:
    if all([n not in u_nodes for n in get_prec(main_graph, scc)]):
        final_sccs.append(scc)

# Printiamo a video il numero di SCC valide del grafo.
print(len(final_sccs))