from generate import *

if __name__ == "__main__":
    
    random.seed(24122022)
    
    GRAPH_NAME = f"sample_test_iter"
    GRAPH_FOLDER = f"{os.getcwd()}/samples/iter_tests/"
    
    for i, CYCLES_TO_CREATE in enumerate(range(100, 25000, 700)):    
        RANDOM_ARCS_TO_ADD = 0.99999
        RANDOM_NODES_TO_ADD = 0.3
        N_ARCS_TO_REMOVE = CYCLES_TO_CREATE*17
        P_SCC_BEING_IN_U = 0.50
        RATIO_U_NODES_TO_TOTAL_NODES = 1/20
        LOWER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*(2/60))
        UPPER_BOUND_LENGTH_CYCLE = int(CYCLES_TO_CREATE*(5/60))
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
        node_of_every_cycle += int(node_of_every_cycle*RANDOM_NODES_TO_ADD)
        main_graph = nx.empty_graph(node_of_every_cycle, create_using=nx.DiGraph)

        print("Unisco i cicli...")
        new_edges = list(map(add_offset_to_edges, cycles))
        del cycles
        new_edges = list(itertools.chain.from_iterable(new_edges))
        
        for _ in range(int(RANDOM_ARCS_TO_ADD*node_of_every_cycle)):
            new_edges.append((random.choice(range(node_of_every_cycle)), random.choice(range(node_of_every_cycle))))

        print(f"Nuovi archi calcolati {len(new_edges)}...")

        main_graph.add_edges_from(new_edges)
        del new_edges
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
        del main_graph
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
        scc_in_u = itertools.chain.from_iterable(random.choices(not_only_one, k=int(P_SCC_BEING_IN_U*len(not_only_one))))
        random_u_nodes = set(random.choices(range(n_nodes),k=int(RATIO_U_NODES_TO_TOTAL_NODES*n_nodes))).union(set(scc_in_u))

        sample_text = f"% {len(order_adj)} {len(final_graph.nodes)}\n"
        sample_text += "\n".join([f"{u} {v}" for u, v in order_adj])
        sample_text += "\n" + "\n".join(list(map(str, random_u_nodes)))
        with open("{}{}{:03}".format(GRAPH_FOLDER, GRAPH_NAME, i), "w") as sample:
            sample.write(sample_text)
