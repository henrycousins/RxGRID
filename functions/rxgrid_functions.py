

import pandas as pd
from collections import Counter
import networkx as nx
import pickle



def get_screenhits(screenname, gene_percentile, rnk_folder, verbose=False):
    """
    Load rnk file
    Return top [gene_percentile] genes from rnk file
    """
    results = pd.read_csv(f'{rnk_folder}/{screenname}.rnk',sep='\t',header=None).sort_values(by=1,ascending=False, kind='mergesort')

    results = list(results[0])
    screenhits = results[0:int(gene_percentile*len(results))]
    if verbose:
        print(screenname)
        print(len(screenhits))
        print(screenhits[0:10] + ['...'])
    return screenhits

def make_full_graph(verbose=False):

    """
    Load drug-chemical interaction network from pickled Nx graph
    Load lists of all genes ("targets") and drugs ("compound") in graph
    """

    G = pickle.load(open('data/G_full.p', 'rb'))
    all_targets = pickle.load(open('data/all_targets.p', 'rb'))
    all_compounds = pickle.load(open('data/all_compounds.p', 'rb'))

    types = nx.get_node_attributes(G, 'type')

    if verbose:
        print('Full graph:')
        print('Number of nodes:',G.number_of_nodes())
        print('Number of edges:',G.number_of_edges())
        for key, value in Counter([types[x] for x in types]).items():
            print(f'Number of type {key}',value)

    return G, types, all_targets, all_compounds


def make_hits_subgraph(G_full, screenhits, remove_singletons=True, verbose=False):
    """
    Remove nodes that are not screenhits
    Remove compounds that have no edges
    Return new graph
    """
    G_hits = G_full.copy() # Start with a copy of the full graph
    G_nodes = [x for x in G_full.nodes()]
    node_types = nx.get_node_attributes(G_hits,'type')
    for i, node in enumerate(G_nodes): # Remove genes that are not screen hits
        t = node_types[node]
        if t == 'target':
            if node not in screenhits:
                G_hits.remove_node(node)
    
    if remove_singletons:
        G_nodes = [x for x in G_hits.nodes()]
        for i, node in enumerate(G_nodes): # Remove drugs that no longer are connected to any genes
            t = node_types[node]
            if t == 'compound':
                if G_hits.degree(node) == 0:
                    G_hits.remove_node(node)

    types = nx.get_node_attributes(G_hits, 'type')

    if verbose:
        print('Hits subgraph:')
        print('Number of nodes:',G_hits.number_of_nodes())
        print('Number of edges:',G_hits.number_of_edges())
        for key, value in Counter([types[x] for x in types]).items():
            print(f'Number of type {key}',value)
        print()

    return G_hits

def compute_compound_normalized_degree_centrality(G_hits, all_compounds):
    """
    Rank nodes by centrality
    Consider only drugs (not genes)
    Sort drugs by centrality score
    Return sorted drugs and scores
    """
    ranking = nx.degree_centrality(G_hits)
    # print(ranking)
    ranking = pd.DataFrame.from_dict(ranking, orient='index')
    ranking = ranking.loc[ranking.index.isin(all_compounds)]
    ranking.rename(columns={0:'score'}, inplace=True)
    ranking.sort_index(inplace=True, ascending=False)
    ranking.sort_values(by='score', ascending=False, inplace=True, kind='stable')
    # print(ranking)
    return ranking


def compute_compound_betweenness_centrality(G_hits, all_compounds):
    """
    Rank nodes by centrality
    Consider only drugs (not genes)
    Sort drugs by centrality score
    Return sorted drugs and scores
    """
    ranking = nx.betweenness_centrality(G_hits)
    ranking = pd.DataFrame.from_dict(ranking, orient='index')
    ranking = ranking.loc[ranking.index.isin(all_compounds)]
    ranking.rename(columns={0:'score'}, inplace=True)
    ranking.sort_values(by='score', ascending=False, inplace=True, kind='mergesort')
    return ranking
