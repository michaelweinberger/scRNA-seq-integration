# A pipeline to integrate scRNA-seq datasets using scVI/scANVI or SAMap
---
Use this pipeline to integrate scRNA-seq datasets using either:
- [scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)<sup>1</sup> and [scANVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scanvi.html)<sup>2</sup>, or
- [SAMap](https://github.com/atarashansky/SAMap?tab=readme-ov-file)<sup>3</sup>

The datasets to be integrated should contain human, mouse or zebrafish data. <br>
Before integration using scVI/scANVI, gene names across datasets will be assimilated with the preference human > mouse > zebrafish (ie. if there is a human dataset, gene names in all other datasets will be converted to human gene names. If there is no human dataset but a mouse and a zebrafish dataset, gene names will be converted to mouse nomenclature). Integration using scVI is not limited to cross-species integration, but can be performed on any metadata variable specified. The output of scVI is used as input to scANVI, which uses cell type labels to optimise data integration. <br>
SAMap specifically works to integrate cross-species datasets, and as a first step generates BLAST databases comparing either whole transcriptomes or coding sequence subsets between species.


## Usage

The pipeline is designed to be run on a high performance server, via the Slurm scheduler. Modules loaded in the scripts that we used to run the analysis are:
R-base/4.3.0, python-base/3.11.3,  blast/2.11.0

Note: The SAMap analysis is run within an Anaconda environment. The following commands have worked for me to set up the environment:
```
$ conda create -n SAMap -c conda-forge python=3.9 numpy=1.23.5 pip pybind11 leidenalg
$ python-igraph texttable
$ conda activate SAMap
$ git clone https://github.com/atarashansky/SAMap.git samap_directory
$ cd samap_directory
$ pip install .
$ pip install --force-reinstall --no-cache-dir holoviews
$ pip install --force-reinstall --no-cache-dir selenium
$ conda install -c conda-forge firefox geckodriver
$ pip install --force-reinstall --no-cache-dir numpy==1.23.5
$ pip install --force-reinstall --no-cache-dir hnswlib==0.7.0
```

To run the pipeline:

1. Clone the repository via:
```
    git clone https://github.com/michaelweinberger/scRNA-seq-integration.git
```
   
2. Adjust the `User defined variables` section of the **1_PARENT_script.sh** script:
### General
- `project`   Name for the project
- `script_dir`   Directory containing scripts copied from https://github.com/michaelweinberger/ scRNA-seq-integration/scripts/
- `out_dir`   Directory containing all output, will be created if non-existent
- `input_file_list`   List of file paths pointing to Scanpy objects containing scRNA-seq data to be integrated. The objects should contain a full matrix of raw expression counts in the .X slot (expression matrix should not be subset to highly variable genes). The .obs metadata in each object should contain an identically named column (e.g. "sample_id" or "cell_type") which will be used in scVI as key for data integration (`scvi_key`) and in SAMap as key to determine the maximum neighbourhood size of each cell (`samap_key`).
- `species_list`   List of species names that the scRNA-seq objects in  `input_file_list` were generated in, one of "human", "mouse" or "zebrafish". The order of `species_list` should match that of `input_file_list`. Each species should be named only once in `species_list`, ie. each input dataset should originate from a different species.

### scVI/scANVI integration
- `scvi_run`   Indicates if integration via scVI and scANVI should be run ("Yes" or "No")
- `scvi_key`   Name of the .obs metadata column in the input Anndata objects that should be used for scVI integration
- `scanvi_key`  Name of the .obs metadata column in the input Anndata objects (column containing cell type labels) that should be used for scanVI integration. Following integration, marker gene identification and benchmarking will be performed on this variable.
- `acc`   Indicates which acceleration should be used ("cpu" or "gpu")
- `n_top_genes`   Indicates the number of highly variable genes to include in the analysis
- `n_layers`   Indicates the number of hidden layers to use for the scVI model
- `n_latent`   Indicates the dimensionality of the latent space
- `scvi_epochs`   Indicates the number of scVI training epochs. If set to "auto", the scVI algorithm will estimate the optimal number of epochs based on dataset size.
- `scanvi_epochs`   Indicates the number of scANVI training epochs
- `min_genes`   Indicates the minimum number of genes detected for a cell to be kept in the datasets
- `max_genes`   Indicates the maximum number of genes detected for a cell to be kept in the datasets
- ` pct_mt_cutoff `   Indicates the maximum percentage of mitochondrial gene counts for a cell to be kept in the datasets
- `min_cells`   Indicates the minimum number of cells in which a gene needs to be detected to be kept in the datasets
- clustering_resolution_scvi `   Indicates the resolution to be used when clustering the integrated dataset

### SAMap integration
- samap_blast   Indicates if transcriptome BLAST databases should be generated ("Yes" or "No"). This step only needs to be performed once for each species.
- samap_run   Indicates if SAMap integration should be run ("Yes" or "No")
- samap_key   Name of the .obs metadata column that should be used for determining maximum neighbourhood size of each cell
- clustering_resolution_samap   Indicates the resolution to be used when clustering the integrated dataset

3. Finally, start the analysis via
```
sbatch 1_PARENT_script.sh
```

## References
1.	Lopez, R., Regier, J., Cole, M.B., Jordan, M.I., and Yosef, N. (2018). Deep generative modeling for single-cell transcriptomics. Nat Methods 15, 1053-1058. 10.1038/s41592-018-0229-2.
2.	Xu, C., Lopez, R., Mehlman, E., Regier, J., Jordan, M.I., and Yosef, N. (2021). Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Mol Syst Biol 17, e9620. 10.15252/msb.20209620.
3.	Tarashansky, A.J., Musser, J.M., Khariton, M., Li, P., Arendt, D., Quake, S.R., and Wang, B. (2021). Mapping single-cell atlases throughout Metazoa unravels cell type evolution. Elife 10. 10.7554/eLife.66747.


