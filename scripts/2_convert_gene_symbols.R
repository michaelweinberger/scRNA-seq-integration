### Run this script to generate a mapping of gene names and IDs between species



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
if (!require("biomaRt")) BiocManager::install("biomaRt")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")



### functions

## function to map gene names and IDs between human, mouse, rat and zebrafish 
# organism_1: organism to map from (hsapiens, mmusculus, rnorvegicus or drerio)
# organism_2: organism to map to (hsapiens, mmusculus, rnorvegicus or drerio)
# out_dir: directory for output file
map_gene_names_species <- function(organism_1, organism_2, out_dir,
                                   out_name_1, out_name_2) {

  # not done due to server errors 
  # read in gtf file 1 (to map from)
  #gtf_df_1 <- as.data.frame(rtracklayer::import(gtf_file_1))

  # extract named genes
  #gtf_df_genes_1 <- gtf_df_1[gtf_df_1$type=="gene",]
  #gtf_df_genes_1 <- gtf_df_genes_1[!is.na(gtf_df_genes_1$gene_name),]
  #genes <- gtf_df_genes_1$gene_id

  #mapping <- getLDS(attributes = c("ensembl_gene_id","external_gene_name"),
  #                  filters = "ensembl_gene_id", values = genes, 
  #                  mart = mart_1,
  #                  attributesL = c("ensembl_gene_id","external_gene_name"), 
  #                  martL = mart_2)


  ## create a mapping of gene id and gene name across species
  if (! file.exists(paste(out_dir, "/jax_gene_mapping_mouse_human_rat_zebrafish.csv", sep=""))) {
    url <- "https://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt"
    destfile <- paste(out_dir, "/HOM_AllOrganism.txt", sep="")

    if (! file.exists(destfile)) download.file(url, destfile)

    mapping_all <- read.table(destfile, header=TRUE, sep="\t", quote="")
    mapping_all <- mapping_all[,c("DB.Class.Key", "Common.Organism.Name",
                                  "Symbol", "EntrezGene.ID")]

    # create dataframe separating entries for separate species into individual columns
    mapping_m <- mapping_all[mapping_all$Common.Organism.Name == "mouse, laboratory",]
    colnames(mapping_m) <- c("DB.Class.Key", "species", "gene_name_mmusculus", "gene_id_mmusculus")
    mapping_m <- mapping_m[,colnames(mapping_m) != "species"]

    mapping_h <- mapping_all[mapping_all$Common.Organism.Name == "human",]
    colnames(mapping_h) <- c("DB.Class.Key", "species", "gene_name_hsapiens", "gene_id_hsapiens")
    mapping_h <- mapping_h[,colnames(mapping_h) != "species"]

    mapping_r <- mapping_all[mapping_all$Common.Organism.Name == "rat",]
    colnames(mapping_r) <- c("DB.Class.Key", "species", "gene_name_rnorvegicus", "gene_id_rnorvegicus")
    mapping_r <- mapping_r[,colnames(mapping_r) != "species"]

    mapping_z <- mapping_all[mapping_all$Common.Organism.Name == "zebrafish",]
    colnames(mapping_z) <- c("DB.Class.Key", "species", "gene_name_drerio", "gene_id_drerio")
    mapping_z <- mapping_z[,colnames(mapping_z) != "species"]

    mapping_m_h <- merge(mapping_m, mapping_h, 
                         by="DB.Class.Key", all=TRUE)
    mapping_m_h_r <- merge(mapping_m_h, mapping_r, 
                           by="DB.Class.Key", all=TRUE)
    mapping <- merge(mapping_m_h_r, mapping_z, 
                     by="DB.Class.Key", all=TRUE)

    write.csv(mapping, file=paste(out_dir, "/jax_gene_mapping_mouse_human_rat_zebrafish.csv", sep=""),
              row.names=FALSE)
  } else {
    mapping <- read.csv(paste(out_dir, "/jax_gene_mapping_mouse_human_rat_zebrafish.csv", sep=""))
  }

  # subset to species of interest
  mapping <- mapping[,c(paste("gene_id_",organism_1,sep=""), 
                        paste("gene_name_",organism_1,sep=""),
                        paste("gene_id_",organism_2,sep=""), 
                        paste("gene_name_",organism_2,sep=""))]

  # simplify column names
  colnames(mapping) <- c("gene_id_1", "gene_name_1",
                         "gene_id_2", "gene_name_2")

  # drop missing entries
  mapping <- mapping[!is.na(mapping$gene_id_1),]
  mapping <- mapping[!is.na(mapping$gene_id_2),]

  # exclude rows without gene names
  rows_drop <- c(rownames(mapping[mapping$gene_name_1=="",]),
                 rownames(mapping[mapping$gene_name_2=="",]))
  mapping <- mapping[! rownames(mapping) %in% rows_drop,]


  ## remove duplicate orthologues (only keep orthologue with highest homology)
  mart_1 <- useMart(biomart = "ensembl", dataset = paste(organism_1, "_gene_ensembl", sep=""),
                    host = "www.ensembl.org/")
  #attributes <- listAttributes(mart_1)

  # ROUND 1:  from organism 1 to organism 2 (one to many)
  genes <- unique(mapping$gene_id_1)

  # get mapping of entrez gene IDs to ensembl gene IDs
  tmp_1 <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                 filters = "entrezgene_id",
                 values = genes, 
                 mart = mart_1)

  # only keep lowest ensembl gene ID if multiple are mapped to single entrez gene ID
  tmp_1 <- tmp_1[order(tmp_1$ensembl_gene_id),]
  tmp_1 <- tmp_1[!duplicated(tmp_1$entrezgene_id),]

  # retrieve gene homology score data
  tmp_2 <- getBM(attributes = c("ensembl_gene_id",
                                paste(organism_2,"_homolog_orthology_confidence",sep=""),
                                paste(organism_2,"_homolog_goc_score",sep=""),
                                paste(organism_2,"_homolog_perc_id_r1",sep="")),
                 filters = "ensembl_gene_id",
                 values = tmp_1$ensembl_gene_id,
                 mart = mart_1)

  # combine entrez gene IDs and gene homology score data
  tmp <- merge(tmp_1, tmp_2, by="ensembl_gene_id", all.x=TRUE)

  # combine mapping and gene homology score data
  mapping_1 <- merge(mapping, tmp, 
                     by.x="gene_id_1",
                     by.y="entrezgene_id", all.x=TRUE)

  # identify duplicate gene entries of organism 1
  mapping_1_no_dup <- as.data.frame(mapping_1 %>% group_by(gene_id_1) %>% filter(n()==1))
  mapping_1_dup <- as.data.frame(mapping_1 %>% group_by(gene_id_1) %>% filter(n()>1))

  # for each duplicated gene entry, check homology data and 
  # choose highest-scoring orthologue in organism 2;
  # if homology data does not allow decision, choose orthologue entry with smaller gene ID;
  # append data to mapping_1_no_dup
  mapping_1_dup <- mapping_1_dup[order(mapping_1_dup$gene_id_2),]
  for (i in unique(mapping_1_dup$gene_id_1)) {
    tmp_df <- mapping_1_dup[mapping_1_dup$gene_id_1 == i,]

    if (length(unique(tmp_df[,paste(organism_2,"_homolog_orthology_confidence",sep="")])) == 1) {
      if (length(unique(tmp_df[,paste(organism_2,"_homolog_goc_score",sep="")])) == 1) {
        if (length(unique(tmp_df[,paste(organism_2,"_homolog_perc_id_r1",sep="")])) == 1) {
          mapping_1_no_dup <- rbind(mapping_1_no_dup, 
                                    tmp_df[1,])
        } else {
          mapping_1_no_dup <- rbind(mapping_1_no_dup, 
                                    tmp_df[tmp_df[,paste(organism_2,"_homolog_perc_id_r1",sep="")] == 
                                                   max(tmp_df[,paste(organism_2,"_homolog_perc_id_r1",sep="")]),][1,])
        }
      } else {
        mapping_1_no_dup <- rbind(mapping_1_no_dup, 
                                  tmp_df[tmp_df[,paste(organism_2,"_homolog_goc_score",sep="")] == 
                                                 max(tmp_df[,paste(organism_2,"_homolog_goc_score",sep="")]),][1,])
      }
    } else {
      mapping_1_no_dup <- rbind(mapping_1_no_dup, 
                                tmp_df[tmp_df[,paste(organism_2,"_homolog_orthology_confidence",sep="")] == 
                                               max(tmp_df[,paste(organism_2,"_homolog_orthology_confidence",sep="")]),][1,])
    }
  }


  ## remove duplicate orthologues (only keep orthologue with highest homology)
  mart_2 <- useMart(biomart = "ensembl", dataset = paste(organism_2, "_gene_ensembl", sep=""),
                    host = "www.ensembl.org/")

  # ROUND 2: from organism 2 to organism 1 (one to many)
  # retrieve gene homology score data
  genes <- unique(mapping_1_no_dup$gene_id_2)

  # get mapping of entrez gene IDs to ensembl gene IDs
  tmp_1 <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                 filters = "entrezgene_id",
                 values = genes, 
                 mart = mart_2)

  # only keep lowest ensembl gene ID if multiple are mapped to single entrez gene ID
  tmp_1 <- tmp_1[order(tmp_1$ensembl_gene_id),]
  tmp_1 <- tmp_1[!duplicated(tmp_1$entrezgene_id),]

  # retrieve gene homology score data
  tmp_2 <- getBM(attributes = c("ensembl_gene_id",
                                paste(organism_1,"_homolog_orthology_confidence",sep=""),
                                paste(organism_1,"_homolog_goc_score",sep=""),
                                paste(organism_1,"_homolog_perc_id_r1",sep="")),
                 filters = "ensembl_gene_id",
                 values = tmp_1$ensembl_gene_id, 
                 mart = mart_2)

  # combine entrez gene IDs and gene homology score data
  tmp <- merge(tmp_1, tmp_2, by="ensembl_gene_id", all.x=TRUE)

  # combine mapping and gene homology score data
  mapping_2 <- merge(mapping_1_no_dup, tmp, 
                     by.x="gene_id_2",
                     by.y="entrezgene_id", all.x=TRUE)

  # identify duplicate gene entries of organism 2
  mapping_2_no_dup <- as.data.frame(mapping_2 %>% group_by(gene_id_2) %>% filter(n()==1))
  mapping_2_dup <- as.data.frame(mapping_2 %>% group_by(gene_id_2) %>% filter(n()>1))

  # for each duplicated gene entry, check homology data and 
  # choose highest-scoring orthologue in organism 1;
  # if homology data does not allow decision, choose orthologue entry with smaller gene ID;
  # append data to mapping_2_no_dup
  mapping_2_dup <- mapping_2_dup[order(mapping_2_dup$gene_id_1),]
  for (i in unique(mapping_2_dup$gene_id_2)) {
    tmp_df <- mapping_2_dup[mapping_2_dup$gene_id_2 == i,]

    if (length(unique(tmp_df[,paste(organism_1,"_homolog_orthology_confidence",sep="")])) == 1) {
      if (length(unique(tmp_df[,paste(organism_1,"_homolog_goc_score",sep="")])) == 1) {
        if (length(unique(tmp_df[,paste(organism_1,"_homolog_perc_id_r1",sep="")])) == 1) {
          mapping_2_no_dup <- rbind(mapping_2_no_dup, 
                                    tmp_df[1,])
        } else {
          mapping_2_no_dup <- rbind(mapping_2_no_dup, 
                                    tmp_df[tmp_df[,paste(organism_1,"_homolog_perc_id_r1",sep="")] == 
                                                   max(tmp_df[,paste(organism_1,"_homolog_perc_id_r1",sep="")]),][1,])
        }
      } else {
        mapping_2_no_dup <- rbind(mapping_2_no_dup, 
                                  tmp_df[tmp_df[,paste(organism_1,"_homolog_goc_score",sep="")] == 
                                                 max(tmp_df[,paste(organism_1,"_homolog_goc_score",sep="")]),][1,])
      }
    } else {
      mapping_2_no_dup <- rbind(mapping_2_no_dup, 
                                tmp_df[tmp_df[,paste(organism_1,"_homolog_orthology_confidence",sep="")] == 
                                               max(tmp_df[,paste(organism_1,"_homolog_orthology_confidence",sep="")]),][1,])
    }
  }


  ## clean up final dataframe + save
  mapping_2_no_dup <- mapping_2_no_dup[,c("gene_id_1", "gene_name_1", "gene_id_2", "gene_name_2")]
  mapping_2_no_dup <- mapping_2_no_dup[order(mapping_2_no_dup$gene_id_1),]
  colnames(mapping_2_no_dup) <- c(paste("gene_id_", out_name_1, sep=""),
                                  paste("gene_name_", out_name_1, sep=""),
                                  paste("gene_id_", out_name_2, sep=""),
                                  paste("gene_name_", out_name_2, sep=""))

  write.csv(mapping_2_no_dup, 
            file=paste(out_dir, "/jax_gene_mapping_", out_name_1, "_", out_name_2, ".csv", sep=""),
            quote=FALSE, row.names=FALSE)

  return(mapping_2_no_dup)
}




### Analysis

organism_list <- c("hsapiens", "mmusculus", "drerio")
out_names_list <- c("human", "mouse", "zebrafish")


# create gene mappings between species
# only run if there are at least 2 organisms in organism_list
if (length(organism_list) > 1) {

  # Initialise list with conditions that remain to be compared
  organism_list_2 <- organism_list

  # do not loop over last organism, because it would not have anything to compare to
  for (i in seq(1,(length(organism_list)-1))) {

    # gradually shrink organism_list_2 to avoid duplicate comparisons
    organism_list_2 <- organism_list_2[organism_list_2 != organism_list[i]]
    print(organism_list_2)

    for (k in seq(1,length(organism_list_2))) {
      if (organism_list[i] != organism_list_2[k]) {
        if (! file.exists(paste(out_dir, "/jax_gene_mapping_", 
                                organism_list[i], "_", organism_list_2[k], ".csv", sep=""))) {
          mapping <- map_gene_names_species(organism_1=organism_list[i],
                                            organism_2=organism_list_2[k],
                                            out_dir=out_dir,
                                            out_name_1=out_names_list[i],
                                            out_name_2=out_names_list_2[k])
        }
      }  
    }
  }
}




