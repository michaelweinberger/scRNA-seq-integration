# -*- coding: utf-8 -*-

"""

Plot SAMap output

Done using Python/11.3

"""




### Modules
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
import argparse

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi_save=450, facecolor='white')
#sc.logging.print_versions()
sns.set_style("ticks")
random.seed(0)



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)
parser.add_argument('-k','--samap_integration_key', help='Name of metadata column to use for samap integration', required=True)
parser.add_argument('-r','--clustering_resolution', help='Clustering resolution of integrated dataset', required=True)

args = parser.parse_args()
out_dir = args.output_dir
out_name = args.output_name
samap_integration_key = args.samap_integration_key
clustering_resolution = float(args.clustering_resolution)



### Functions

def obs_bar(adata, column_1, column_2, out_dir, out_name):
    
    # count cell type labels by leiden cluster
    plot_df = adata.obs.groupby([column_1, column_2]).agg({column_2: 'size'})
    plot_df = plot_df.rename(columns={column_2:'count'})
    plot_df = plot_df.reset_index()
    
    # generate dataframe with percentage counts for each entry in column_2 in separate column
    plot_df_1 = plot_df.pivot(index=column_1, columns=column_2, values='count')   
    plot_df_1 = plot_df_1.div(plot_df_1.sum(axis=1), axis=0)
    
    # generate stacked bar plot
    values = plot_df_1.to_dict(orient='list')
    labels = plot_df_1.index
    bottom = np.zeros(len(labels))
    width = 0.8
    #fig_width = len(labels) * 0.7
    
    fig, ax = plt.subplots(figsize=(14,9), frameon=False)
    
    for key, value in values.items():
        p = ax.bar(labels, value, width, label=key, bottom=bottom)
        bottom += value
        #ax.bar_label(p, label_type='center')
        
    ax.set_axisbelow(True)
    ax.set_title("")
    ax.set_xlabel(column_1, fontsize=13)
    ax.set_ylabel("Relative count", fontsize=13)
    ax.tick_params(axis='x', which='major', labelsize=13, labelrotation=30)
    ax.tick_params(axis='y', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title=column_2, 
              title_fontsize=12, fontsize=13, frameon=False)
    
    ax.figure.savefig(f"{out_dir}/{out_name}_{column_1}_vs_{column_2}_stacked_bar.pdf",
                      dpi=300)
    
    plt.close()
    return(plot_df_1)




def samap_plot(samap_obj, leiden_res, out_dir, out_name):
    
    """
    Function to plot UMAP plots using SAMap output
    
    Input:
    samap_obj: Anndata object generated via SAMap
    leiden_res: Leiden clustering resolution
    out_dir: Directory to write Anndata object and plots to
    out_name: Name of dataset, will be used in plot and output file names
    
    """
    
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # set plot file suffix
    sc.settings.plot_suffix = f"_{out_name}"
    sc.settings.figdir = out_dir

    # plot stacked bar charts
    sc.tl.leiden(samap_obj, resolution=leiden_res)
    
    #print(samap_obj)
    
    plot1 = obs_bar(adata=samap_obj, column_1="leiden", column_2="dataset", 
                   out_dir=out_dir, out_name=f"{out_name}_samap")
    
    plot2 = obs_bar(adata=samap_obj, column_1=samap_integration_key, column_2="dataset", 
                   out_dir=out_dir, out_name=f"{out_name}_samap")
    
    # drop unnamed columns from metadata
    samap_obj.obs = samap_obj.obs.loc[:, ~samap_obj.obs.columns.str.contains('^Unnamed')]
 
    # plot umap with metadata overlay
    for col in samap_obj.obs.columns:
        if col != 'barcode':
            if col != 'barcode_tmp':
                sc.pl.umap(samap_obj, color=[col], 
                           legend_loc='right margin', title='', frameon=False, ncols=1, 
                           save=f"_{col}_samap.png")
    
    return(samap_obj)





### Analysis
    
samap_obj = sc.read_h5ad(f"{out_dir}/{out_name}_samap.h5ad")
    
samap_plot(samap_obj=samap_obj, leiden_res=clustering_resolution, out_dir=out_dir, out_name=out_name)


