#!/bin/env python3

import pygraphviz as pgv
import argparse

### EDGES TO FIGURE ###
# Simple script to parse output haplotype edges and draw network diagram from the output of the AnVir pipeline.
# NOTE: The output file is read and parsed twice. The first parse just grabs the metadata required to filter
# out nodes worth including in the diagram. The second read/parse reads all the edge data to construct the
# network.

# AnVir Haplotype edges file configuration constants
HEADER_SIZE = 2

# INPUT COLUMN VALUES
# ID - int - the unique node id, note: not necessarily consistant between differnt files/runs.
ID = 0

# COUNT - int - number of observations of the haplotypes.
COUNT = 1

# BITSET - str - a list of the variants that define the haplotype e.g. { 1 3 17 }
BITSET = 2

# STATUS - str/enum - the existance of some nodes is inferred rather than observed directly, this status indicates whether
# the node haplotype has been observed directlyi (key = orig) or not.
STATUS = 3

# N_CHILDREN - int - the number of outgoing edges from the node.
N_CHILDREN = 4

# N_ANCESTORS - int - the number of inbound edges towards node.
N_ANCESTORS = 5

# FINAL COLUMNS(3) REPEAT - each repeated 3 column block represents a single edge.
# The first set of edges are children/decendents, the second set of edges are the ancestors.
# There is an "X" to mark the transition.
EDGES_START = 6
EDGES_INTERVAL = 3

EDGE_ID = 0
EDGE_BITSET = 1
EDGE_POSTERIOR = 2 # theoretically this is the posterior likelihood value of the edge, but don't trust them at all.

# PARSER 1 - retrieves metadata for nodes from input file.
def parse_haplo_node_row(row_elements):
    """
    Each set of row elements represents a single haplotype node.
    """
    return({
        "ID" : int(row_elements[ID]), 
        "count" : int(row_elements[COUNT]),
        "bitset" : row_elements[BITSET],
        "status" : row_elements[STATUS],
        "n_children" : int(row_elements[N_CHILDREN]),
        "n_ancestors" : int(row_elements[N_ANCESTORS])
    })
    
def read_nodes_from_file(file_name):
    """
    Parse input haplotype edges file to get metadata.
    Returns list of dicts representing each node (keys are column names)
    """
    with open(file_name, "r") as fh:
        lines = [line for line in fh.readlines()][HEADER_SIZE:]
        haplo_nodes = [parse_haplo_node_row(line.strip().split("\t")) for line in lines]
        return(haplo_nodes)

# PARSER 2 - retrieves the all the edges between nodes from input file.
def parse_children_and_ancestors(edge_elements):
    """
    Returns an EDGES DICT representing the edges connected to a single haplotype node.
    Contains two list (CHILDREN, ANCESTORS) of dicts (bitset, id) which represent the nodes at the other end of
    the each edge.
    """
    edges_at_node = {
        "CHILDREN" : [],
        "ANCESTORS" : []
    }
    
    state = "CHILDREN"
    while(len(edge_elements) > 0):
        if(edge_elements[0] == "X"):
            # FOUND "X" marker - add edges to Ancestors rather than children.
            state = "ANCESTORS"
            edge_elements = edge_elements[1:]
            continue
        
        edges_at_node[state].append({ "bitset": edge_elements[EDGE_BITSET],
                                      "id" : edge_elements[EDGE_ID] })
        edge_elements = edge_elements[EDGES_INTERVAL:]
    return(edges_at_node)

def read_edges_from_data(target_haplos, file_name):
    """
    Parses the edge infomation from the input file for only the target nodes.
    target_haplos = set of bitsets representing the haplotypes to be retrieved from the input file.
    Returns dict mapping nodes to their edges (NODE_BITSET -> EDGES DICT)
    """
    nodes_to_edges = {}
    with open(file_name, "r") as fh:
        lines = [line for line in fh.readlines()][HEADER_SIZE:]
        for i, line in enumerate(lines):
            row_elements = line.strip().split("\t")
            if row_elements[BITSET] in target_haplos:
                nodes_to_edges[row_elements[BITSET]] = parse_children_and_ancestors(row_elements[EDGES_START:])
    return(nodes_to_edges)

# Construct graphviz graph necessary to draw output figure.
def build_haplo_graphviz(all_nodes, node_edges):
    """
    Returns graphviz datastructure representing haplotype network.
    Basic formatting infomation is included e.g. node size.
    """
    BASE_NODE_SIZE = 0.02
    FREQ_SIZE_SCALE = 1.2
    SIZE_CUTOFF = 0.8 # nodes above this adjusted size with have their bitsets labeled in figure.

    # Max observed frequency of haplotype - used to find upper limit of node size.
    max_count = max([x["count"] for x in all_nodes.values()])

    G = pgv.AGraph(strict=False, directed=True)
    G.graph_attr["outputorder"] = "edgesfirst"
    G.node_attr["shape"] = "circle"
    G.node_attr["penwidth"] = 1.0
    G.edge_attr["penwidth"] = 0.7
    G.edge_attr["arrowsize"] = 0.2
    G.edge_attr["color"] = "gray"
    G.node_attr["fixedsize"] = True
    
    for node_bitset, node_metadata in all_nodes.items():
        if((FREQ_SIZE_SCALE * (node_metadata["count"]/float(max_count))) > SIZE_CUTOFF):
            G.add_node(node_bitset, height = BASE_NODE_SIZE + (FREQ_SIZE_SCALE * (node_metadata["count"]/float(max_count))))
        else:
            G.add_node(node_bitset, height = BASE_NODE_SIZE + (FREQ_SIZE_SCALE * (node_metadata["count"]/float(max_count))), label="")
    for node_bitset, edges in node_edges.items():
        for child in edges["CHILDREN"]:
            if child["bitset"] in all_nodes:
                G.add_edge(node_bitset, child["bitset"])
    return(G)

if __name__ =="__main__":
    # Parse Command Line Arguments
    parser = argparse.ArgumentParser(description="Draws a basic network diagram from AnVir output file.")
    parser.add_argument("input_file", help="AnVir haplotype edges file (often .xls)")
    parser.add_argument("output_file", help="given file extension will determine file type e.g. .png or .svg")

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    print(f"INPUT: {input_file}")

    # Read/parse for metadata.
    all_haplo_nodes = read_nodes_from_file(input_file)

    # Filter nodes to those to be included in network diagram.
    # Removes all nodes not observed directly or that have no edges associated with them.
    target_nodes = {node["bitset"]: node for node in
                   filter(lambda node: (node["status"] == "orig") and
                          ((node["n_children"] > 0) or (node["n_ancestors"] > 0)),
                          all_haplo_nodes)}

    print(f"Constructing network containing {len(target_nodes)} haplotype node.")

    # Read/paarse for edge infomation.
    nodes_to_edges = read_edges_from_data(set(target_nodes.keys()), input_file)

    # Build network
    network = build_haplo_graphviz(target_nodes, nodes_to_edges)

    print(f"OUTPUT: {output_file}")
    network.layout(prog="neato")
    network.draw(output_file)
