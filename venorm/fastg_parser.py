#!/usr/bin/env python3

import networkx
import matplotlib
from Bio import SeqIO


def fetch_subgraph_contigs(contigs_parent_dir, paths, i=1, visualise=False):
    '''
    Fetch any contigs with connectivity to the longest assembly contig by parsing FASTG output.
    Assumes that the longest contig is the first contig, which it is for SPAdes.
    Canonicalise forward and reverse complement nodes
    Handle cases where longest contig is unconnected
    Returns subgraph of contigs and its number of nodes
    '''
    graph = networkx.Graph()
    longest_contig = None
    longest_contig_unconnected = None
    headers_cnames = {}
    with open(contigs_parent_dir + '/contigs.fastg', 'r') as contigs_fastg_file:
        for record in SeqIO.parse(contigs_fastg_file, 'fasta'): # treat fastg as fasta
            node_name, node_neighbors = None, None
            canonicalised_header = record.id[:-1].replace("'","")
            if ':' in canonicalised_header: # is node connected?
                node_name, node_neighbors = canonicalised_header.split(':')
                if longest_contig is None:
                   longest_contig = node_name
                node_neighbors = node_neighbors.split(',')
                for node_neighbor in node_neighbors:
                    if (node_name, node_neighbor) not in graph.edges():
                        graph.add_edge(node_name, node_neighbor)
            else:
                node_name = canonicalised_header
            if longest_contig is None:
                longest_contig = node_name
                graph.add_node(node_name)

    subgraph = graph.subgraph(networkx.node_connected_component(graph, longest_contig))
    subgraph_nodes = subgraph.nodes()
    subgraph_node_lens = [int(node.split('length_')[1].split('_cov')[0]) for node in subgraph_nodes]
    subgraph_nodes_lens = dict(zip(subgraph_nodes, subgraph_node_lens))

    with open(contigs_parent_dir + '/contigs.fasta', 'r') as contigs_file:
        with open(paths['o'] + '/remap/subgraph_contigs.fasta', 'w') as subgraph_contigs_file:
            for record in SeqIO.parse(contigs_file, 'fasta'):
                if record.id in subgraph_nodes:
                    SeqIO.write(record, subgraph_contigs_file, 'fasta')

    subgraph_node_labels = [str(item) for item in subgraph_node_lens]
    subgraph_nodes_labels = dict(zip(subgraph_nodes, subgraph_node_labels))
    positions = networkx.spring_layout(subgraph)

    if visualise:
        networkx.draw(subgraph, pos=positions, node_size=subgraph_node_lens, with_labels=False)
        networkx.draw_networkx_labels(subgraph, pos=positions, labels=subgraph_nodes_labels)
        matplotlib.pyplot.show()

    return subgraph, len(subgraph_nodes)

# contigs_parent_dir = '/Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp/run_1443791863_phe_sample8/asm/1.Sample8.HIV-generic.1-0..norm_k25c2.asm_k21,33,55,77.rg0/'
# fetch_subgraph_contigs(contigs_parent_dir, paths)
