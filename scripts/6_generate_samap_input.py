# -*- coding: utf-8 -*-
"""


"""

### Modules
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import os
import argparse


## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)

args = parser.parse_args()
out_dir = args.output_dir



### Analysis

species_list = open(f"{out_dir}/species.txt").read().splitlines()
dataset_list = open(f"{out_dir}/datasets.txt").read().splitlines()

for dataset,species in zip(dataset_list, species_list):

    out_dir_tmp = os.path.dirname(dataset)
    out_name_tmp = os.path.splitext(os.path.basename(dataset))[0]
    
    # adjust object for SAMap analysis
    if not os.path.isfile(f"{out_dir_tmp}/{out_name_tmp}_samap_input.h5ad"):
        
        adata = sc.read_h5ad(dataset)
        adata.obs['dataset'] = out_name_tmp
        adata.obs['barcode'] = adata.obs.index
        adata.obs.index = adata.obs['barcode'] + "_" + adata.obs['dataset']
    
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
    
        adata.write(f"{out_dir_tmp}/{out_name_tmp}_samap_input.h5ad", compression='gzip')
