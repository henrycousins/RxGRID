import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib.patches import Rectangle
import argparse
import time

from functions.rxgrid_functions import get_screenhits, make_full_graph, make_hits_subgraph, compute_compound_normalized_degree_centrality, compute_compound_betweenness_centrality


if __name__ == '__main__':

    print('Running RxGRID...')

    # LOAD COMMAND-LINE ARGUMENTS
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_hit_threshold', type=float, required=False, default=0.05,
                        help='Significance cutoff for genes (default=0.05)') # Upper proportion of genes considered hits
    parser.add_argument('--drug_hit_threshold', type=float, required=False, default=0.01,
                        help='Significance cutoff for drugs (default=0.01)') # Upper proportion of drugs considered hits
    parser.add_argument('--centrality_metric', type=str, required=False, default='ndc',
                        help='Primary centrality metric for ranking drugs (degree ["ndc"] or betweenness ["btw"]; default="ndc"') # Centrality metric to use
    parser.add_argument('--mutual_hit_freq', required=False, default=None,
                        help='Fraction of screens in which drug must be hit to meet overall significance (inputs >=1 are interpreted as counts, rather than fractions; default = 1/3)') # Number of screens in which drug should be hit to meet aggregate hit criteria
    parser.add_argument('--drug_exclusions_path', required=False, default='data/drug_exclusions.csv',
                        help='CSV file containing user-specified list of drugs to exclude from downstream analysis (default="data/drug_exclusions.csv")') # CSV file containing list of drugs to exclude from downstream analysis
    parser.add_argument('--rnk_folder', type=str, required=False, default='screens',
                        help='Name of folder containing the screens, as rnk files (default="screens")') # Folder containing screen files
    parser.add_argument('--output_folder', type=str, required=False, default='results',
                        help='Name of folder to which results should be saved (default="results")') # Existing folder to which results should be saved
    parser.add_argument('--save_plots', type=bool, required=False, default=False,
                        help='Whether to save descriptive plots (default=False)') # Whether to save plots
    args = parser.parse_args()

    gene_percentile = args.gene_hit_threshold
    drug_percentile = args.drug_hit_threshold
    metric = args.centrality_metric
    rnk_folder = args.rnk_folder
    save_to = args.output_folder
    save_plots = args.save_plots
    drug_exclusions_path = args.drug_exclusions_path

    timestamp = str(int(time.time()))

    # Screen names are automatically derived from file names
    # for x in os.listdir(rnk_folder):
    #     print()
    #     print(x)
    #     assert x[-4:] == '.rnk', f'Found filenames in screen directory missing required ".rnk" extension'
    screennames = [x[0:-4] for x in os.listdir(rnk_folder) if x[-4:]=='.rnk']
    print(f'Identified {len(screennames)} rnk files: {[x+".rnk" for x in screennames]}')

    # Check that centrality metric is supported
    assert metric in ['ndc','btw'], 'Currently supported centrality metrics are "ndc" [normalized degree centrality] and "btw" [betweenness centrality]'

    # Required number of screen hits will be forced as integer
    required_number_of_screen_hits = args.mutual_hit_freq 
    if required_number_of_screen_hits == None:
        required_number_of_screen_hits = int(np.round(len(screennames)/3)) # Default
    elif required_number_of_screen_hits < 1:
        required_number_of_screen_hits = int(np.round(len(screennames)*required_number_of_screen_hits)) # If input is a fraction, treat it as such
    else:
        required_number_of_screen_hits = int(np.round(required_number_of_screen_hits)) # Input will be left as is if it is an integer greater than one

    # User-defined list of drugs to exclude from analysis (such as narcotics, dietary supplements, etc.)
    if drug_exclusions_path == None:
        drugs_to_exclude = []
    else:
        drugs_to_exclude = list(pd.read_csv(drug_exclusions_path, header=None)[0])



    # GENERATE FULL GRAPH
    G_full, types, all_targets, all_compounds = make_full_graph()


    # GENERATE SCREEN-SPECIFIC SUBGRAPHS

    graph_stats = pd.DataFrame()
    drug_dict = {}
    drug_dict_info = {}
    drug_screen_hits = pd.DataFrame()
    screen_drug_values = {}

    for screenname in screennames: # For each screen

        print(f'Processing {screenname}...')

        screenhits = get_screenhits(screenname, gene_percentile, rnk_folder, verbose=False) # Get gene hits
        G_hits = make_hits_subgraph(G_full, screenhits, verbose=False, remove_singletons=True) # Make subgraph based on screen hits

        # Calculate descriptive statistics
        num_genes = len([x for x in G_hits.nodes if x in all_targets]) # number of genes in subgraph
        num_drugs = len([x for x in G_hits.nodes if x in all_compounds]) # number of drugs in subgraph
        num_edges = len([x for x in G_hits.edges]) # number of edges in subgraph
        mean_degree = np.mean([G_hits.degree(x) for x in G_hits.nodes]) # mean degree of nodes in subgraph
        mean_degree_genes = np.mean([G_hits.degree(x) for x in G_hits.nodes if x in all_targets]) # mean degree of genes in subgraph
        mean_degree_drugs = np.mean([G_hits.degree(x) for x in G_hits.nodes if x in all_compounds]) # mean degree of drugs in subgraph
        density = num_edges / (num_genes * num_drugs) # density

        graph_stats.loc[screenname,'num_genes'] = num_genes
        graph_stats.loc[screenname,'num_drugs'] = num_drugs
        graph_stats.loc[screenname,'num_edges'] = num_edges
        graph_stats.loc[screenname,'mean_degree'] = mean_degree
        graph_stats.loc[screenname,'mean_degree_genes'] = mean_degree_genes
        graph_stats.loc[screenname,'mean_degree_drugs'] = mean_degree_drugs
        graph_stats.loc[screenname,'density'] = density

        # This procedure ensures stable ranking:
        betweenness_ranking = compute_compound_betweenness_centrality(G_hits, all_compounds)
        ndc_ranking = compute_compound_normalized_degree_centrality(G_hits, all_compounds)
        agg_ranking = pd.merge(ndc_ranking, betweenness_ranking, left_index=True, right_index=True,
                            suffixes=('_ndc','_btw'))
        agg_ranking.sort_index(ascending=True, inplace=True)
        agg_ranking.sort_values(by='score_btw', ascending=False, inplace=True, kind='mergesort')
        if metric == 'ndc':
            agg_ranking.sort_values(by='score_ndc', ascending=False, inplace=True, kind='mergesort')
        if metric == 'btw':
            agg_ranking.sort_values(by='score_btw', ascending=False, inplace=True, kind='mergesort')
        
        # For this screen, drug hits are defined as follows:
        top_drugs = list(agg_ranking.iloc[0:int(np.ceil(len(all_compounds)*drug_percentile))].index)
        drug_dict[screenname] = top_drugs # Save the collection of drug hits
        screen_drug_values[screenname] = agg_ranking.copy() # Save the full ranking



    # SAVE SUBGRAPH DESCRIPTIONS

    graph_stats.to_csv(f'{save_to}/{timestamp}_graph_descriptions.csv')


    # COMBINE SCREEN-SPECIFIC HITS

    drug_hits_combined = []
    for screenname in drug_dict.keys():
        drug_hits_combined += drug_dict[screenname]
    drug_hits_count = Counter(drug_hits_combined)

    collective_hits = []
    for drug in drug_hits_count.keys():
        if drug_hits_count[drug] >= required_number_of_screen_hits:
            collective_hits.append(drug)

    final_drugs = [x for x in collective_hits if x not in drugs_to_exclude]

    print(f'Identified {len(final_drugs)} significant drugs:')
    print(sorted(final_drugs))

    # Save results

    print('Saving results...')

    graph_stats = graph_stats.transpose()
    graph_stats['mean'] = graph_stats.apply(lambda x: np.nanmean(x[screennames]), axis=1)
    graph_stats['standev'] = graph_stats.apply(lambda x: np.nanstd(x[screennames]), axis=1)

    graph_stats.to_csv(f'{save_to}/{timestamp}_graph_descriptions.csv')


    # Plotting


    ####

    save_figure3_plots = False
    plt.rcParams['font.size'] = str(16)

    ####


    df_results = pd.DataFrame(columns=[
        'Screenwise hit count',
        'Centrality score mean',
        'Centrality score standard deviation',
        'Percentile rank mean',
        'Percentile rank standard deviation',
        ] + [f'{screenname} is_hit' for screenname in screennames] + [f'{screenname} percentile_rank' for screenname in screennames] + [f'{screenname} centrality_score' for screenname in screennames],

        )


    heatframe = pd.DataFrame()
    boxdata = []
    for screenname in screennames:

        screen_ranking = screen_drug_values[screenname]
        screen_ranking['rank'] = pd.DataFrame([x for x in range(screen_ranking.shape[0])], index=screen_ranking.index) # Drug rank raw
        screen_ranking['percentile_rank'] = screen_ranking['rank'].apply(lambda x: (len(all_compounds) - x)/len(all_compounds) * 100) # Drug rank, converted to percentile
        screen_ranking['centrality_score'] = screen_ranking[f'score_{metric}'].apply(lambda x: (x/max(screen_ranking[f'score_{metric}']))) # Drug centrality score

        for drug in sorted(final_drugs):
            if drug in screen_ranking.index:
                drug_percentile_rank = screen_ranking['percentile_rank'].loc[drug]
                drug_centrality_score = screen_ranking['centrality_score'].loc[drug]
            else:
                drug_percentile_rank = np.nan
                drug_centrality_score = 0
            boxdata.append([drug, drug_centrality_score])
            heatframe.loc[screenname, drug] = drug_percentile_rank

            if drug in drug_dict[screenname]:
                drug_is_hit = 1
            else:
                drug_is_hit = 0

            df_results.loc[drug, f'{screenname} is_hit'] = drug_is_hit
            df_results.loc[drug, f'{screenname} percentile_rank'] = drug_percentile_rank
            df_results.loc[drug, f'{screenname} centrality_score'] = drug_centrality_score


    df_results['Screenwise hit count'] = df_results.apply(lambda x: sum(x[[f'{screenname} is_hit' for screenname in screennames]]), axis=1)
    df_results['Percentile rank mean'] = df_results.apply(lambda x: np.nanmean(x[[f'{screenname} percentile_rank' for screenname in screennames]]), axis=1)
    df_results['Percentile rank standard deviation'] = df_results.apply(lambda x: np.nanstd(x[[f'{screenname} percentile_rank' for screenname in screennames]]), axis=1)
    df_results['Centrality score mean'] = df_results.apply(lambda x: np.nanmean(x[[f'{screenname} centrality_score' for screenname in screennames]]), axis=1)
    df_results['Centrality score standard deviation'] = df_results.apply(lambda x: np.nanstd(x[[f'{screenname} centrality_score' for screenname in screennames]]), axis=1)

    df_results.sort_values(by='Screenwise hit count', ascending=False, inplace=True)

    df_results.to_csv(f'{save_to}/{timestamp}_summary.csv')


    if save_plots:

        print('Creating plots...')

        boxdata = pd.DataFrame(boxdata, columns=['drug','centralityscore'])

        fig, ax = plt.subplots()
        fig.set_figheight(3)
        fig.set_figwidth(12)
        fig.set_figwidth(len(final_drugs)*0.5)

        b = sns.boxplot(x='drug', y='centralityscore', data=boxdata, 
                    color='lightgray',
                    ax=ax,
                    )

        ax.set_ylabel('Centrality score')

        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
        ax.set_xlabel('')

        for box in b.artists:
            box.set_edgecolor('k')
        plt.setp(b.lines, color='k')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        plt.savefig(f'{save_to}/{timestamp}_centralityplot.png', dpi=1000, bbox_inches='tight')


        fig,ax = plt.subplots()
        fig.set_figheight(3)
        # fig.set_figwidth(15)
        fig.set_figwidth(len(final_drugs)*0.5)

        h = sns.heatmap(heatframe,
                        linecolor = 'k',
                        linewidths = 1,
                        cmap='Blues',
                        cbar_kws={
                            'label': 'Percentile rank',
                            'pad': 0.01,
                            },
                        ax=ax)

        ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

        for i,screenname in enumerate(heatframe.index):
            for j,drugname in enumerate(heatframe.columns):
                if drugname in drug_dict[screenname]:
                    ax.add_patch(Rectangle((j, i), 1, 1, fill=False, edgecolor='red', lw=1, hatch='//'))

        plt.savefig(f'{save_to}/{timestamp}_screenwiseplot.png', dpi=1000, bbox_inches='tight')

print('RxGRID complete.')