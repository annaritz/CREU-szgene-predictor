import networkx as nx


def main():

    GIANT = nx.read_weighted_edgelist("brain_top.txt")
    GIANT_nodes = GIANT.nodes()
    print(type(GIANT_nodes))
    print("Number of GIANT nodes: ",len(GIANT_nodes))

    HUGO_nodes = read_edge_file("hgnc.txt")
    print("Number of HUGO nodes: ", len(HUGO_nodes))

    compare_nodes(GIANT_nodes,HUGO_nodes)


def read_edge_file(edge_file):
    nodes = set()
    with open(edge_file, "r") as file:
        for line in file: 
            node = line.strip().split("\t") #strips whitespace characters
            nodes.add(node[1]) #adds the entrez gene ID to the set of nodes
    return nodes


def compare_nodes(Gset, Hset):
    intersection = Gset.intersection(Hset) #creates a new set of the intersection (common elements) of the sets
    print("Number of common nodes:", len(intersection))

    Gdiff = Gset.difference(Hset) #creates new set of the nodes in GIANT but not in HUGO
    print("Number of nodes in GIANT but not in HUGO:",len(Gdiff))

    Hdiff = Hset.difference(Gset) #creates new set of the nodes in HUGO but not in GIANT
    print("Number of nodes in HUGO but not in GIANT:",len(Hdiff))
    return


if __name__ == "__main__":
    main()