#!/bin/bash


#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G

#SBATCH --partition=gpu
#SBATCH --gpus=2
 


#########################################################
### This script:
# integrates scRNA-seq datasets across species
#########################################################



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

###################################################################################################

# specify project name
project="Kuppe_Wang"

# specify directory containing parent and child scripts
script_dir="/ceph/project/tsslab/mweinber/tmp/scripts_scvi"

# specify output directory
# if directory does not exist, it will be generated
out_dir="/ceph/project/tsslab/mweinber/tmp/out_scvi"

# supply list of file paths to non-clustered Scanpy objects containing count matrix of all genes
input_file_list=("/ceph/project/tsslab/mweinber/2023_datasets/datasets_scRNAseq/Kuppe_et_al_10.1038s41586-022-05060-x/scanpy/Kuppe_scRNAseq_no_doublets_annotated_scVI.h5ad"
"/ceph/project/tsslab/mweinber/2023_datasets/datasets_scRNAseq/Wang_et_al_10.1016j.celrep.2020.108472/scanpy/Wang_scRNAseq_no_doublets_annotated_scVI.h5ad"
)

# specify list of species that datasets have been generated in ("human", "mouse" or "zebrafish"),
# order needs to match datasets supplied in input_file_list
species_list=("human" "mouse")



##############   scVI + scANVI integration   ##############

# indicate if scVI and scANVI should be run ("Yes" or "No")
scvi_run="Yes"

# indicate name of metadata column that should be used for scVI integration
scvi_key="sample_id"

# indicate name of metadata column (containing cell type labels) that should be used for scanVI integration
# Following integration, marker gene identification and benchmarking will be performed on this variable
scanvi_key="cell_type"

# indicate which acceleration is used ("cpu" or "gpu")
acc="gpu"

# indicate the number of highly variable genes to include in analysis
n_top_genes=5000

# indicate the number of hidden layers to use for the scVI model
n_layers=3

# indicate the dimensionality of the latent space
n_latent=20

# indicate the number of scVI training epochs (if set to "auto", the algorithm will estimate the optimal number of epochs based on dataset size)
scvi_epochs="auto"

# indicate the number of scANVI training epochs
scanvi_epochs=50

# indicate the minimum number of genes a cell needs to express to be kept in the datasets
min_genes=400

# indicate the maximum number of genes a cell needs to express to be kept in the datasets
max_genes=5000

# indicate the minimum number of cells a gene needs to be expressed in to be kept in the datasets
min_cells=3

# indicate the maximum fraction of mitochondrial reads allowed for a cell to be kept in the datasets
pct_mt_cutoff=15

# indicate cell clustering resolution of integrated dataset
clustering_resolution_scvi=0.4



##############   SAMap integration   ##############

# indicate if transcriptome blast databases should be generated ("Yes" or "No")
# only needs to be done once for each species
samap_blast="No"


# indicate if SAMap integration should be run ("Yes" or "No")
samap_run="No"

# indicate name of metadata column that should be used for determining maximum neighbourhood size of each cell
samap_key="cell_type"

# indicate cell clustering resolution of integrated dataset
clustering_resolution_samap=0.4


###################################################################################################





### Analysis

module load R-cbrg
module load python-cbrg
module load blast/2.11.0

cd ${script_dir}



## create output directory if it does not exist
[ ! -d "$out_dir" ] && mkdir -p "$out_dir"



## download genome gtf files
len=${#species_list[@]}

for (( k=0; k<${len}; k++ ))
do

	# specify additional genome information
	if test "${species_list[k]}" = "human" ; then
		species_latin="homo_sapiens"
		species_latin_2="Homo_sapiens"
		species_latin_short="hsapiens"
		genome="GRCh38"
		genome_ucsc="hg38"
		chrom_number=22 #X,Y chromosomes excluded
	elif test "${species_list[k]}" = "mouse" ; then
		species_latin="mus_musculus"
		species_latin_2="Mus_musculus"
		species_latin_short="mmusculus"
		genome="GRCm39"
		genome_ucsc="mm39"
		chrom_number=19 #X,Y chromosomes excluded
	elif test "${species_list[k]}" = "zebrafish" ; then
		species_latin="danio_rerio"
		species_latin_2="Danio_rerio"
		species_latin_short="drerio"
		genome="GRCz11"
		genome_ucsc="danRer11"
		chrom_number=25
	fi

	# save species and genome name arrays, to be used in R script
	if test $k = 0;then
		echo ${input_file_list[k]} > ${out_dir}/datasets.txt
		echo ${species_list[k]} > ${out_dir}/species.txt
		echo $genome > ${out_dir}/genomes.txt
		echo $genome_ucsc > ${out_dir}/genomes_ucsc.txt
		echo $species_latin_short > ${out_dir}/species_latin_short.txt
		echo $chrom_number > ${out_dir}/chrom_number.txt
	else
		echo ${input_file_list[k]} >> ${out_dir}/datasets.txt
		echo ${species_list[k]} >> ${out_dir}/species.txt
		echo $genome >> ${out_dir}/genomes.txt
		echo $genome_ucsc >> ${out_dir}/genomes_ucsc.txt
		echo $species_latin_short >> ${out_dir}/species_latin_short.txt
		echo $chrom_number >> ${out_dir}/chrom_number.txt
	fi
done



## integrate data via scVI + scANVI
if test $scvi_run = "Yes"; then

	# generate gene mappings between species
	out_dir_gene_mappings=${out_dir}/gene_mappings
	[ ! -d "$out_dir_gene_mappings" ] && mkdir -p "$out_dir_gene_mappings"

	R -f ${script_dir}/2_convert_gene_symbols.R --args out_dir=${out_dir_gene_mappings}

	# run scVI integration
	out_dir_scvi=${out_dir}/scVI
	[ ! -d "$out_dir_scvi" ] && mkdir -p "$out_dir_scvi"

	python ${script_dir}/3_scVI.py -o $out_dir_scvi -n $project -g $out_dir_gene_mappings -s ${out_dir}/species.txt -d ${out_dir}/datasets.txt \
		-k $scvi_key -kk $scanvi_key -a $acc -tg $n_top_genes -l $n_layers -lt $n_latent -e $scvi_epochs -ee $scanvi_epochs \
		-mng $min_genes -mxg $max_genes -mnc $min_cells -mt $pct_mt_cutoff -r $clustering_resolution_scvi
fi



## integrate data via SAMap
out_dir_tr=${out_dir}/SAMap/transcriptomes
[ ! -d "$out_dir_tr" ] && mkdir -p "$out_dir_tr"

out_dir_cds=${out_dir}/SAMap/coding_sequences
[ ! -d "$out_dir_cds" ] && mkdir -p "$out_dir_cds"


# generate transcriptome blast databases for species
if test $samap_blast = "Yes"; then

	# generate transcriptome and coding sequence fasta files
	R -f ${script_dir}/4_make_transcriptome_fasta.R --args in_dir=${out_dir} out_dir_tr=${out_dir_tr} out_dir_cds=${out_dir_cds}

	# generate blast comparisons between fasta files
	len=${#species_list[@]}

	for (( i=0; i<${len}; i++ ))
	do
		species_1=${species_list[i]}

		if test "${species_1}" = "human" ; then
			n1="hs"
		elif test "${species_1}" = "mouse" ; then
			n1="mm"
		elif test "${species_1}" = "zebrafish" ; then
			n1="dr"
		fi

		for (( k=${i}; k<${len}; k++ ))
		do
			species_2=${species_list[k]}

			if test "${species_2}" = "human" ; then
				n2="hs"
			elif test "${species_2}" = "mouse" ; then
				n2="mm"
			elif test "${species_2}" = "zebrafish" ; then
				n2="dr"
			fi

			# create blast maps from transcriptome fasta files
			if [ ! -f "${out_dir_tr}/maps/${n1}${n2}/${n1}_to_${n2}.txt" ]; then
				cd ${out_dir_tr}
				# map_genes.sh downloaded from https://github.com/atarashansky/SAMap/blob/main/
				source ${script_dir}/5_map_genes.sh --tr1 ${out_dir_tr}/${species_1}_transcriptome.fa \
             	         	 	--t1 "nucl" --n1 ${n1} --tr2 ${out_dir_tr}/${species_2}_transcriptome.fa \
             		         	--t2 "nucl" --n2 ${n2}
			fi

			# create blast maps from coding sequence fasta files
			if [ ! -f "${out_dir_cds}/maps/${n1}${n2}/${n1}_to_${n2}.txt" ]; then
				cd ${out_dir_cds}
				# map_genes.sh downloaded from https://github.com/atarashansky/SAMap/blob/main/
				source ${script_dir}/5_map_genes.sh --tr1 ${out_dir_cds}/${species_1}_cds.fa \
                       			--t1 "nucl" --n1 ${n1} --tr2 ${out_dir_cds}/${species_2}_cds.fa \
                       			--t2 "nucl" --n2 ${n2}
			fi
		done
	done
fi


# run SAMap
if test $samap_run = "Yes"; then

	# generate SAMap input objects
	python ${script_dir}/6_generate_samap_input.py -o $out_dir

	# run SAMap
	module unload python-cbrg
	conda run -n SAMap python ${script_dir}/7_samap.py -o $out_dir -n $project -k $samap_key -r $clustering_resolution_samap
	module load python-cbrg

	# plot UMAPs using SAMap output
	python ${script_dir}/8_samap_plotting.py -o ${out_dir}/SAMap/out/transcriptome/${project} -n $project -k $samap_key -r $clustering_resolution_samap
	python ${script_dir}/8_samap_plotting.py -o ${out_dir}/SAMap/out/coding_sequences/${project} -n $project -k $samap_key -r $clustering_resolution_samap
fi



echo All done!

