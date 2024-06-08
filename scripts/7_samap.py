# -*- coding: utf-8 -*-

"""

Integrate scRNA-seq datasets via SAMap

Done using Python/3.9

"""


#pip install samap --user --no-dependencies
#pip install sam-algorithm==1.0.2 --user --no-dependencies
#pip install hnswlib==0.7.0 --user --no-dependencies


### Modules
import numpy as np
import pandas as pd
import re as re
import matplotlib.pyplot as plt
import seaborn as sns
import os
import random
import samap as samap
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import holoviews as hv
import argparse

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

def samap_run(adata_dict, transcriptome_map_dir, alignment_keys, clustering_dict, 
              resolutions, out_dir, out_name):
    
    """
    Function to integrate scRNA-seq data via SAMap
    
    Input:
    adata_dict: Dictionary specifying filepaths of input Anndata objects, 
                keys: species IDs determined in the BLAST step
    transcriptome_map_dir: Directory containing Blast-compared transcriptome files
    alignment_keys: Dictionary of metadata column names in Anndata objects to use for SAMap integration,
                    keys: species IDs determined in the BLAST step
    clustering_dict: Dictionary of clustering resolutions to use for clustering of integrated dataset,
                     keys: species IDs determined in the BLAST step
    resolutions: None or dictionary of leiden clustering resolutions in Anndata objects,
                 keys: species IDs determined in the BLAST step
    out_dir: Directory to write output to
    out_name: Name of dataset, will be used in plot and output file names
    
    """
    
    hv.extension('bokeh', logo=False)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # generate SAMap object
    print(f"Generating SAMap object using {transcriptome_map_dir}")    
    sm = SAMAP(adata_dict, f_maps=f"{transcriptome_map_dir}/", keys=alignment_keys, resolutions=resolutions)
    
    sm.run(pairwise=True)
    samap_obj = sm.samap # SAM object with input datasets stitched together
    
    # compute alignment scores
    D,MappingTable = get_mapping_scores(sm, alignment_keys, n_top = 0)
    
    # plot alignment scores
    fig = sankey_plot(MappingTable, align_thr=0.05)
    hv.save(fig, f"{out_dir}/{out_name}_samap_alignment_sankey.png", fmt="png")
    
    fig = chord_plot(MappingTable, align_thr=0.05)
    hv.save(fig, f"{out_dir}/{out_name}_samap_alignment_chord.png", fmt="png")
    
    # plot cell clustering
    fig = sm.scatter(ss=clustering_dict)
    fig.figure.savefig(f"{out_dir}/{out_name}_samap_clustering.png", dpi=300)    
    
    # find enriched gene pairs across alignment key categories
    gpf = GenePairFinder(sm, keys=alignment_keys)
    gene_pairs = gpf.find_all(align_thr=0.2)
    gene_pairs.to_csv(f"{out_dir}/{out_name}_enriched_gene_pairs.csv")
    

    ## recombine columns that have been split across species
    # drop unnamed columns from metadata
    samap_obj.adata.obs = samap_obj.adata.obs.loc[:, ~samap_obj.adata.obs.columns.str.contains('^Unnamed')]
        
    # convert categorical metadata into strings
    for col in samap_obj.adata.obs.columns:
        if samap_obj.adata.obs.dtypes[col] == "category":
            samap_obj.adata.obs[col] = samap_obj.adata.obs[col].astype(str)
        
    # replace string value with NaN
    samap_obj.adata.obs.replace('unassigned', np.nan, inplace=True)
    
    # split pre-pended species ID from original column name
    split_cols = list()
    
    for i in range(0, len(samap_obj.adata.obs.columns)):
        if len(samap_obj.adata.obs.columns[i].split('_')) > 1:
            if samap_obj.adata.obs.columns[i].split('_', 1)[0] in ["hs", "mm", "dr"]:
                split_cols.append(samap_obj.adata.obs.columns[i].split('_', 1)[1])
        
    # identify duplicate column names
    dup_cols = {x for x in split_cols if split_cols.count(x) > 1}
    
    # recombine duplicate columns
    for col in dup_cols:
        
        # get columns across species
        col_list = samap_obj.adata.obs.filter(regex=(f".*{col}$")).columns.to_list()
        
        # create new column in which to recombine columns
        samap_obj.adata.obs[col] = samap_obj.adata.obs[col_list[0]]
        
        # fill NA in new column with corresponding entries from other species column(s)
        for i in range(1, len(col_list)):
            samap_obj.adata.obs[col] = samap_obj.adata.obs[col].fillna(samap_obj.adata.obs[col_list[i]])
            
        # drop original columns
        samap_obj.adata.obs.drop(col_list, axis=1, inplace=True)

    samap_obj.adata.obs.fillna("unassigned", inplace=True)

    # adjust alignment key column created by SAMap
    samap_obj.adata.obs.columns = samap_obj.adata.obs.columns.str.replace(";", "_")
    old_name = samap_obj.adata.obs.filter(regex=(f".*_mapping_scores$")).columns.to_list()[0]
    samap_obj.adata.obs[old_name] = samap_obj.adata.obs[old_name].str.split('_', n=1).str[1]     
    
    new_name = alignment_keys[next(iter(alignment_keys))]
    samap_obj.adata.obs[new_name] = samap_obj.adata.obs[old_name]
    samap_obj.adata.obs.drop([old_name], axis=1, inplace=True)

    # save data files
    samap.utils.save_samap(sm, f"{out_dir}/{out_name}_samap_obj")
    samap_obj.adata.write(f"{out_dir}/{out_name}_samap.h5ad", compression='gzip')
    
    return(samap_obj)





### Analysis

transcriptome_map_dir = f"{out_dir}/SAMap/transcriptomes/maps"
cds_map_dir = f"{out_dir}/SAMap/coding_sequences/maps"

species_list = open(f"{out_dir}/species.txt").read().splitlines()
dataset_list = open(f"{out_dir}/datasets.txt").read().splitlines()

dataset_dict = dict()

for dataset,species in zip(dataset_list, species_list):

    out_dir_tmp = os.path.dirname(dataset)
    out_name_tmp = os.path.splitext(os.path.basename(dataset))[0]
    
    # create dictionary of input files
    if species == "human":
        dataset_dict["hs"] = f"{out_dir_tmp}/{out_name_tmp}_samap_input.h5ad"
    elif species == "mouse":
        dataset_dict["mm"] = f"{out_dir_tmp}/{out_name_tmp}_samap_input.h5ad"
    elif species == "zebrafish":
        dataset_dict["dr"] = f"{out_dir_tmp}/{out_name_tmp}_samap_input.h5ad"


# create dictionary of metadata variables on which to align datasets
alignment_keys = dataset_dict
alignment_keys = dict.fromkeys(alignment_keys, samap_integration_key)

# create dictionary of cell clustering resolutions
clustering_dict = dataset_dict
clustering_dict = dict.fromkeys(clustering_dict, clustering_resolution)


# run SAMap using whole gene BLAST databases
samap_run(adata_dict=dataset_dict, 
          transcriptome_map_dir=transcriptome_map_dir, 
          alignment_keys=alignment_keys, 
          clustering_dict=clustering_dict, 
          resolutions=None,
          out_dir=f"{out_dir}/SAMap/out/transcriptome/{out_name}", 
          out_name=out_name
)

# run SAMap using coding sequence BLAST databases
samap_run(adata_dict=dataset_dict, 
          transcriptome_map_dir=cds_map_dir, 
          alignment_keys=alignment_keys, 
          clustering_dict=clustering_dict,
          resolutions=None,
          out_dir=f"{out_dir}/SAMap/out/coding_sequences/{out_name}", 
          out_name=out_name
)

