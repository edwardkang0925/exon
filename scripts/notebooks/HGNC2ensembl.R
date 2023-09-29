library(data.table)
library(dplyr)
library(tidyr)  # For the gather function used later
library(GENESIS)

gtf.file.path <- "data/gencode.v38.annotation.gtf"

genes <- fread(gtf.file.path)
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

# [optional] focus, for example, only on entries of type "gene",
# which will drastically reduce the file size
genes <- genes[type == "gene"]

# the problem is the attributes column that tends to be a collection
# of the bits of information you're actually interested in
# in order to pull out just the information I want based on the
# tag name, e.g. "gene_id", I have the following function:
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

# this is how to, for example, extract the values for the attributes of interest (here: "gene_id")
genes$gene_id <- unlist(lapply(genes$attributes, extract_attributes, "gene_id"))
genes$gene_name <- unlist(lapply(genes$attributes, extract_attributes, "gene_name"))

# create HGNC -> Ensembl mapping
map <- genes[,c('gene_id', 'gene_name')]

# save it
saveRDS(map, "data/HGNC2ensembl.RDS")

# checking for preprocessForConditional_all_updating.R
ensembl_id = 'ENSG00000076770.15'
# given gene_id, convert to hgnc using map then check for bins containing the gene under gene_name metadata column.
mapped_hgnc_name = map[map$gene_id == ensembl_id,]$gene_name
mapped_hgnc_name
subsetBins <- bins[sapply(elementMetadata(bins)$gene_name, function(x) mapped_hgnc_name %in% strsplit(as.character(x), split=",")[[1]])]
# create a vector of bin name
bin_ids <- paste(elementMetadata(subsetBins)$gene, elementMetadata(subsetBins)$bin_id, sep=".")
bin_ids

# Given gene name in ENSG format, return all binID within the gene.
geneToBins <- function(gene_ensembl_id, bins, map) {
  # Convert ensembl_id to hgnc_name using the map
  mapped_hgnc_name <- map[map$gene_id == gene_ensembl_id, ]$gene_name

  # Subset the bins based on the mapped_hgnc_name
  subsetBins <- bins[sapply(elementMetadata(bins)$gene_name, function(x)
    mapped_hgnc_name %in% strsplit(as.character(x), split=",")[[1]])
  ]

  # Create a vector of bin names
  bin_ids <- paste(elementMetadata(subsetBins)$gene, elementMetadata(subsetBins)$bin_id, sep=".")

  # append the gene too.
  bin_ids <- c(bin_ids, gene_ensembl_id)

  return(bin_ids)
}

geneToBins(ensembl_id, bins, map)

