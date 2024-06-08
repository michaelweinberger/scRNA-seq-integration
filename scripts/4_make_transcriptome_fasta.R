### Run this script to generate a mapping of genes between species



### user defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)
args

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}





### packages

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")
if (!require("BSgenome")) BiocManager::install("BSgenome")
if (!require("BSgenome.Hsapiens.UCSC.hg38")) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
if (!require("BSgenome.Mmusculus.UCSC.mm39")) BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
if (!require("BSgenome.Drerio.UCSC.danRer11")) BiocManager::install("BSgenome.Drerio.UCSC.danRer11")
if (!require("gprofiler2")) BiocManager::install("gprofiler2")
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene")) BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene", lib="/ceph/home/m/mweinber/R_packages")
library("TxDb.Mmusculus.UCSC.mm39.knownGene", lib="/ceph/home/m/mweinber/R_packages")





### functions

## function to generate nucleotide transcriptome fasta file annealing exon sequences for each gene
make_transcriptome_fasta <- function(organism, chrom_num, dataset, genome, out_dir, out_name) {
  
  if (organism == "hsapiens") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (organism == "mmusculus") {
    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  } else {
    txdb <- makeTxDbFromBiomart(dataset=dataset)

    # convert chromosome names to UCSC format
    newSeqNames <- paste('chr', seqlevels(txdb), sep = '')
    names(newSeqNames) <- seqlevels(txdb)
    txdb <- renameSeqlevels(txdb, newSeqNames)
  }

  # subset to main chromosomes
  seqlevels(txdb) <- paste("chr", seq(1, chrom_num), sep="")

  # generate GRanges list object, listing exons per gene
  tmp <- exonsBy(txdb, by="gene")

  # combine overlapping exons for individual genes
  tmp1 <- reduce(tmp)

  # get DNA sequences
  fasta <- extractTranscriptSeqs(genome, tmp1)

  # convert gene IDs to gene symbols
  ids2names <- gconvert(query = fasta@ranges@NAMES, organism = organism, 
                        target="ENSG", numeric_ns="ENTREZGENE_ACC", 
                        mthreshold = Inf, filter_na = TRUE)
  ids2names <- ids2names[ids2names$name!="None",]
  ids2names <- ids2names[!duplicated(ids2names$input),]

  # find gene symbols of gene IDs in sequence object
  fasta_names <- as.data.frame(fasta@ranges@NAMES)
  colnames(fasta_names) <- "gene_id"
  fasta_id_names <- fasta_names %>% left_join(ids2names[,c("input","name")], 
                                              by=join_by(gene_id==input))

  # check that order of gene ID to gene symbol file matches that of gene IDs in sequence object
  print(paste("New names are identical to old names: ", identical(fasta_id_names$gene_id,
                                                                  fasta@ranges@NAMES), sep=""))

  # replace gene IDs in sequence object with gene symbols
  fasta@ranges@NAMES <- fasta_id_names$name

  # drop genes without names
  fasta <- fasta[!is.na(fasta@ranges@NAMES),]

  # save as fasta file
  writeXStringSet(fasta, filepath=paste(out_dir, "/", out_name, "_transcriptome.fa", sep=""),
                  format="fasta")
  return(fasta)             
}



## function to generate nucleotide coding sequence fasta file 
make_cds_fasta <- function(organism, chrom_num, dataset, genome, out_dir, out_name) {
  
  if (organism == "hsapiens") {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (organism == "mmusculus") {
    txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
  } else {
    txdb <- makeTxDbFromBiomart(dataset=dataset)

    # convert chromosome names to UCSC format
    newSeqNames <- paste('chr', seqlevels(txdb), sep = '')
    names(newSeqNames) <- seqlevels(txdb)
    txdb <- renameSeqlevels(txdb, newSeqNames)
  }

  # subset to main chromosomes
  seqlevels(txdb) <- paste("chr", seq(1, chrom_num), sep="")

  # generate GRanges list object, listing coding sequences per gene
  tmp <- cdsBy(txdb, by="gene")

  # combine overlapping exons for individual genes
  tmp1 <- reduce(tmp)

  # get DNA sequences
  fasta <- extractTranscriptSeqs(genome, tmp1)

  # convert gene IDs to gene symbols
  ids2names <- gconvert(query = fasta@ranges@NAMES, organism = organism, 
                        target="ENSG", numeric_ns="ENTREZGENE_ACC", 
                        mthreshold = Inf, filter_na = TRUE)
  ids2names <- ids2names[ids2names$name!="None",]
  ids2names <- ids2names[!duplicated(ids2names$input),]

  # find gene symbols of gene IDs in sequence object
  fasta_names <- as.data.frame(fasta@ranges@NAMES)
  colnames(fasta_names) <- "gene_id"
  fasta_id_names <- fasta_names %>% left_join(ids2names[,c("input","name")], 
                                              by=join_by(gene_id==input))

  # check that order of gene ID to gene symbol file matches that of gene IDs in sequence object
  print(paste("New names are identical to old names: ", identical(fasta_id_names$gene_id,
                                                                  fasta@ranges@NAMES), sep=""))

  # replace gene IDs in sequence object with gene symbols
  fasta@ranges@NAMES <- fasta_id_names$name

  # drop genes without names
  fasta <- fasta[!is.na(fasta@ranges@NAMES),]

  # save as fasta file
  writeXStringSet(fasta, filepath=paste(out_dir, "/", out_name, "_cds.fa", sep=""),
                  format="fasta")
  return(fasta)             
}





### Analysis
species <- unlist(read.table(paste(in_dir, "/species.txt", sep=""), sep=" "))
genomes <- unlist(read.table(paste(in_dir, "/genomes_ucsc.txt", sep=""), sep=" "))
organisms <- unlist(read.table(paste(in_dir, "/species_latin_short.txt", sep=""), sep=" "))
chrom_nums <- unlist(read.table(paste(in_dir, "/chrom_number.txt", sep=""), sep=" "))


for (i in seq(1,length(species))) {

  if (species[i] == "human") {
    dataset <- "hsapiens_gene_ensembl"
    genome <- BSgenome.Hsapiens.UCSC.hg38
  } else if (species[i] == "mouse") {
    dataset <- "mmusculus_gene_ensembl"
    genome <- BSgenome.Mmusculus.UCSC.mm39
  } else if (species[i] == "zebrafish") {
    dataset <- "drerio_gene_ensembl"
    genome <- BSgenome.Drerio.UCSC.danRer11
  }
  
  out_name <- species[i]
  
  print(paste("Generating transcriptome fasta file for ", out_name, sep="")) 
  
  if (!file.exists(paste(out_dir_tr, "/", out_name, "_transcriptome.fa", sep=""))) {
    make_transcriptome_fasta(organism=as.list(organisms[i])[[1]],
                             chrom_num=as.list(chrom_nums[i])[[1]],
                             dataset=dataset, genome=genome,
                             out_dir=out_dir_tr, out_name=out_name)
  }

  print(paste("Generating coding sequence fasta file for ", out_name, sep="")) 
  
  if (!file.exists(paste(out_dir_cds, "/", out_name, "_cds.fa", sep=""))) {
    make_cds_fasta(organism=as.list(organisms[i])[[1]],
                            chrom_num=as.list(chrom_nums[i])[[1]],
                            dataset=dataset, genome=genome,
                            out_dir=out_dir_cds, out_name=out_name)
  }
}




