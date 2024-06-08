# -*- coding: utf-8 -*-

"""

Integrate scanpy object datasets with scVI and scanVI

Done using Python/3.11.3

"""



### Modules
# !pip install --quiet scvi-colab
# from scvi_colab import install
# install()

#import matplotlib
#from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import random
import scvi
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker
import argparse

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi_save=300, facecolor='white')
#sc.logging.print_versions() 
random.seed(0)

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
torch.set_float32_matmul_precision("high")



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)
parser.add_argument('-g','--gene_mappings_dir', help='File path to a directory containing gene mappings between species', required=True)
parser.add_argument('-s','--species_list', help='File path to a blank-separated text file listing species analysed', required=True)
parser.add_argument('-d','--dataset_list', help='File path to a blank-separated text file listing datasets to be integrated', required=True)
parser.add_argument('-k','--scvi_integration_key', help='Name of metadata column to use for scvi integration', required=True)
parser.add_argument('-kk','--scanvi_integration_key', help='Name of cell type labels metadata column to use for scanvi integration', required=True)
parser.add_argument('-a','--acceleration', help='Whether to use CPU or GPU acceleration', required=True)
parser.add_argument('-tg','--top_genes', help='Number of highly variable genes to use', required=True)
parser.add_argument('-l','--layers', help='Number of hidden layers to use for the scVI model', required=True)
parser.add_argument('-lt','--latent', help='Dimension of the latent space', required=True)
parser.add_argument('-e','--scvi_epochs', help='Number of scVI training epochs', required=True)
parser.add_argument('-ee','--scanvi_epochs', help='Number of scanVI training epochs', required=True)
parser.add_argument('-mng','--min_genes', help='Minimum number of genes expressed in a valid cell', required=True)
parser.add_argument('-mxg','--max_genes', help='Maximum number of genes expressed in a valid cell', required=True)
parser.add_argument('-mnc','--min_cells', help='Minimum number of cells for a valid gene to be expressed in', required=True)
parser.add_argument('-mt','--pct_mt', help='Maximum fraction of mitochondrial reads in a valid cell', required=True)
parser.add_argument('-r','--clustering_resolution', help='Clustering resolution of integrated dataset', required=True)


args = parser.parse_args()
out_dir = args.output_dir
out_name = args.output_name
gene_mappings_dir = args.gene_mappings_dir
species_list_file = args.species_list
dataset_list_file = args.dataset_list
scvi_integration_key = args.scvi_integration_key
scanvi_integration_key = args.scanvi_integration_key
acceleration = args.acceleration
top_genes = int(args.top_genes)
layers = int(args.layers)
latent = int(args.latent)
scanvi_epochs = int(args.scanvi_epochs)
min_genes = int(args.min_genes)
max_genes = int(args.max_genes)
min_cells = int(args.min_cells)
pct_mt = int(args.pct_mt)
clustering_resolution = float(args.clustering_resolution)

if args.scvi_epochs == "auto":
    scvi_epochs = None
else:
    scvi_epochs = int(args.scvi_epochs)




### Functions

def scvi_integrate(adata, 
                   leiden_res, batch_key, labels_key, acc, 
                   out_dir, out_name, n_top_genes=5000, 
                   n_layers=3, n_latent=20, scvi_epochs=None, scanvi_epochs=50,
                   min_genes=400, min_cells=3, 
                   n_genes_cutoff=5000, pct_mt_cutoff=15):
    
    '''
    
    
    
    '''
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    # set plot file suffix
    sc.settings.plot_suffix = f"_{out_name}"
    sc.settings.figdir = out_dir
       
    # filter data
    sc.pl.highest_expr_genes(adata, n_top=20, save='.png')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # quality control
    # analyse mitochondrial gene content
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, save='_qc_violin.png')
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_qc_scatter_1.png')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_qc_scatter_2.png')
    
    adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
    adata = adata[adata.obs.pct_counts_mt < pct_mt_cutoff, :]
    
    # store raw counts in "counts" layer
    adata.layers['counts'] = adata.X.copy()
    
    # normalise data
    sc.pp.normalize_total(adata, target_sum=1e5)
    sc.pp.log1p(adata)
    #print(adata.X.sum(axis=1))
    #print(adata.layers['counts'].sum(axis=1))    
    
    adata.raw = adata  # keep full dimension safe
    
    # subset to highly variable genes
    sc.pp.highly_variable_genes(adata, 
                                n_top_genes=n_top_genes, span=0.5,
                                batch_key=batch_key, subset=True)
    sc.pl.highly_variable_genes(adata, save='.png')
    
    # eliminate unwanted variability
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # scale gene expression to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)
    
     
    # integrate data with scVI
    #print(adata.layers['counts'].sum(axis=1))
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_layers=n_layers, n_latent=n_latent, gene_likelihood="nb")
    model.train(max_epochs=scvi_epochs, accelerator=acc)
    
    # get scVI corrected counts for HVGs
    adata.layers['scvi_counts'] = model.get_normalized_expression(library_size="latent")
    
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(adata, resolution=leiden_res)
    
    # embed using mde
    SCVI_MDE_KEY = "X_scVI_MDE"
    adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])
    
    # embed using umap
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=True, edge_width_scale=2, node_size_scale=2, save='_scVI.pdf')
    
    sc.tl.umap(adata, init_pos='paga')  
    
    # drop unnamed columns from metadata
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains('^Unnamed')]

    # plot with metadata overlay
    for col in adata.obs.columns:
        if col != 'barcode':
            if col != 'barcode_tmp':
                
                sc.pl.embedding(adata, basis=SCVI_MDE_KEY, color=[col], frameon=False, ncols=1,
                                save=f"_{col}_scVI.png")
                
                sc.pl.umap(adata, color=[col], 
                           legend_loc='right margin', title='', frameon=False, ncols=1, 
                           save=f"_{col}_scVI.png")
    
    # integrate data with scanVI (using cell labels)
    scanvi_model = scvi.model.SCANVI.from_scvi_model(model, adata=adata, labels_key=labels_key,
                                                     unlabeled_category="Unknown")
    
    scanvi_model.train(max_epochs=scanvi_epochs, n_samples_per_label=None, accelerator="gpu")
    
    SCANVI_LATENT_KEY = "X_scANVI"
    adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
    
    # embed using mde
    SCANVI_MDE_KEY = "X_scANVI_MDE"
    adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY])
    
    # embed using umap
    sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
    sc.tl.leiden(adata, resolution=leiden_res)    
    #sc.tl.paga(adata)
    #sc.pl.paga(adata, plot=True, edge_width_scale=2, node_size_scale=2, save='_scVI.pdf')
    
    sc.tl.umap(adata)
    
    # plot with metadata overlay
    for col in adata.obs.columns:
        if col != 'barcode':
            if col != 'barcode_tmp':
            
                sc.pl.embedding(adata, basis=SCANVI_MDE_KEY, color=[col], frameon=False, ncols=1,
                                save=f"_{col}_scanVI.png")
                
                sc.pl.umap(adata, color=[col], 
                           legend_loc='right margin', title='', frameon=False, ncols=1, 
                           save=f"_{col}_scanVI.png")
    
    # identify marker genes
    sc.tl.rank_genes_groups(adata, labels_key, method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='.pdf')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    markers = pd.DataFrame({group + '_' + key[:1]: result[key][group] 
                            for group in groups for key in ['names', 'pvals']})
    markers.to_csv(f"{out_dir}/{out_name}_markers.csv")
    
   
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}_scVI_scanVI.h5ad", compression='gzip')
    
    
    # compute integration metrics
    bm = Benchmarker(adata, batch_key=batch_key, label_key=labels_key,
                     embedding_obsm_keys=["X_pca", SCVI_LATENT_KEY, SCANVI_LATENT_KEY],
                     n_jobs=-1)
    bm.benchmark()
    
    df = bm.get_results(min_max_scale=False).transpose()
    df.to_csv(f"{out_dir}/{out_name}_scvi_scanvi_benchmarking.csv")
    
    bm.plot_results_table(min_max_scale=False, show=False, save_dir=out_dir)
    
    print(adata)
    return(adata)




def assimilate_gene_names(adata, species, target_species, gene_mappings_file):
    
    gene_mappings = pd.read_csv(gene_mappings_file)
    adata.var['features'] = adata.var_names

    tmp_df = adata.var.merge(gene_mappings, 
                             left_on="features", right_on=f"gene_name_{species}", 
                             how="left", suffixes=(None,"_1"), sort=False)
    tmp_df.index = tmp_df['features']
    check_order = adata.var.index.equals(tmp_df.index)
    print(f".var order after adding metadata matches: {check_order}")       
    adata.var = tmp_df
    adata = adata[:,adata.var[f"gene_name_{target_species}"].notna()]
    adata.var_names = adata.var[f"gene_name_{target_species}"]

    adata.var = adata.var[['features']]
    adata.var.index.name = None
    
    print(adata)
    return(adata)






### Analysis
    
species_list = open(species_list_file).read().splitlines()
dataset_list = open(dataset_list_file).read().splitlines()


# identify species to convert gene annotations to
if "human" in species_list:
    target_species = "human"
else:
    if "mouse" in species_list:
        target_species = "mouse"
    else:
        target_species = "zebrafish"


# adjust gene annotations across datasets + concatenate datasets
if not os.path.isfile(f"{out_dir}/{out_name}_scvi_input.h5ad"):
    counter = 0
    for dataset,species in zip(dataset_list, species_list):
        
        print(f"Processing {dataset}")
        adata_tmp = sc.read_h5ad(dataset)
        adata_tmp.obs['dataset'] = os.path.splitext(os.path.basename(dataset))[0]
        adata_tmp.obs['barcode'] = adata_tmp.obs.index
        adata_tmp.obs.index = adata_tmp.obs['barcode'] + "_" + adata_tmp.obs['dataset']
    
        # adjust gene annotation
        if species != target_species:
            print(f"Converting gene names from {species} to {target_species}")
            adata_tmp = assimilate_gene_names(adata=adata_tmp, species=species, target_species=target_species,
                                              gene_mappings_file=f"{gene_mappings_dir}/jax_gene_mapping_{target_species}_{species}.csv")
            #print(adata_tmp.var)
        
        # concatenate datasets
        if counter == 0:
            adata = adata_tmp
            counter = 1
        else: 
            adata = ad.concat([adata, adata_tmp], join='outer', merge='first')
            adata = adata[:, adata.var['features'].notna()]
            adata.var_names_make_unique()

    adata.write(f"{out_dir}/{out_name}_scvi_input.h5ad", compression='gzip')
else: 
    adata = sc.read_h5ad(f"{out_dir}/{out_name}_scvi_input.h5ad")

print(adata)


# integrate with scVI + scanVI
scvi_integrate(adata=adata, batch_key=scvi_integration_key, labels_key=scanvi_integration_key, acc=acceleration,
               leiden_res=clustering_resolution, out_dir=out_dir, out_name=out_name,
               n_top_genes=top_genes, n_layers=layers, n_latent=latent, 
               scvi_epochs=scvi_epochs, scanvi_epochs=scanvi_epochs,
               min_genes=min_genes, min_cells=min_cells, 
               n_genes_cutoff=max_genes, pct_mt_cutoff=pct_mt)


