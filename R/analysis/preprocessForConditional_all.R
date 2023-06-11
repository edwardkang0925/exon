.libPaths(c("/home/edwardkang/R/x86_64-conda-linux-gnu-library/4.2", .libPaths()))

library(dplyr)
library(tibble)
library(stats)

### Functions
mapBinToGene <- function(bin_id) {
  # Remove the part of the string after the last period
  gene_id <- sub("\\.[^.]*$", "", bin_id)
  
  # If PAR_Y is present, remove it
  gene_id <- sub("_PAR_Y", "", gene_id)
  
  return(gene_id)
}

# calculate GIF
inflation <- function(p) {
  chisq <- qchisq(1 - p, df = 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

#######

# this script does basic EPA thing on exon and gene level summary statistics. write logs in a file specified by $logPath
# Also, output merged slice of exon and gene level residuals for conditional analysis.
# this version is preparing inputs for conditional analysis on all bins (if not all are within a gene)
# Let's log how many bins are not contained in a gene. (expecting low count)

# Obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the first argument to trait
trait <- args[1]
numSlice = 100
logPath <- paste0("outputs/log_preprocessForConditional_allBins_",trait,".txt")
FDR = 0.01
write(paste0("----------",trait,"------------"), logPath, append = TRUE)
# output path
output.conditional.genesis.dir <- file.path("outputs", "conditionalGenesisInput_all", trait)

# create output dirs if not existing. 
if(!dir.exists(output.conditional.genesis.dir)){
  dir.create(output.conditional.genesis.dir, recursive=TRUE)
}

# necessary files (TWAS results)
exon_level_summary_file_path <- file.path("outputs", "combined", trait, paste0("combinedGenesis_", trait, ".csv"))
exon_level_residual_file_path <- file.path("outputs", "afterRegression", trait, "combined", "residualExonExpression.RDS")
sigGene_level_summary_file_path <- file.path("data", "geneLvMMAPresults", 
                                             paste0("parsed_output_twas_input_visit_1_adjusted_", trait,
                                                    "_z_scores_p_values_corrected_significant_genes_bonferroni.csv"))
gene_level_residual_file_path <- file.path("data", "MMAP_INPUT_FINAL_ROUND", paste0("twas_input_visit_1_adjusted_",trait,".csv"))
# Read summary data of a trait SLOW
exon.pvals.df <- read.csv(exon_level_summary_file_path)
gene.pvals.df <- read.csv(sigGene_level_summary_file_path)
exon.resid.df <- readRDS(exon_level_residual_file_path)
gene.resid.df <- read.csv(gene_level_residual_file_path)

# bin-level significance threshold
exon.level.sig.threshold <- FDR / dim(exon.resid.df)[2]

# collect significant genes from gene-level TWAS 
# TODO: ENSG00000185291.12 is duplicated for lnTG. probably, PAR_Y thing. check this in exon -> exon data retains PAR_Y
# For now, using unique(), print number of sig bins 
allGenes <- unique(colnames(gene.resid.df)[grepl("ENSG", colnames(gene.resid.df))])
sigGenes <- unique(gene.pvals.df$Genes)
sigExon.df <- exon.pvals.df[exon.pvals.df$pval_bacon < exon.level.sig.threshold,]
write(paste(trait, "has", length(sigGenes), "significant genes"), logPath, append = TRUE)
write(paste(trait, "has", dim(sigExon.df)[1], "significant bins. bacon_p < ", FDR, "/ #bins"), logPath, append = TRUE)


# print GIF 
GIF_nominal <- inflation(exon.pvals.df$pval)
write(paste(trait, "bin-level GIF (nominal) =", GIF_nominal), logPath, append = TRUE)

GIF_bacon <- inflation(exon.pvals.df$pval_bacon)
write(paste(trait, "bin-level GIF (bacon) =", GIF_bacon), logPath, append = TRUE)

# Logging purpose Check exon -> gene coverage & sig-gene coverage
sigExon.df <- sigExon.df %>%
  mutate(gene = mapBinToGene(exon))
genes.containing.sigExons <- unique(sigExon.df$gene)
write(paste(trait, "#unique genes containing sig bin:", length(genes.containing.sigExons)), logPath, append=TRUE)
sigGenesCoveredBysigBins <- intersect(genes.containing.sigExons, sigGenes) # just for logging
write(paste(trait, "#unique significant genes containing sig bin:", length(sigGenesCoveredBysigBins)), logPath, append=TRUE)

# identify all genes containing any bin. (preparing input for conditional analysis on all bins within a gene)
genesWithBin <- unique(lapply(colnames(exon.resid.df), mapBinToGene))
allGenesCoveredByanyBins <- intersect(allGenes, genesWithBin)
# setdiff(genesWithBin, allGenesCoveredByanyBins) # moment I realized what I was thinking as "gene expression" was "transcript expression" .. 


# Subset exon.resid.df to only include these columns
binToGene <- sapply(colnames(exon.resid.df), mapBinToGene)
isMapped.vector <- binToGene %in% allGenesCoveredByanyBins
exon.resid.subset.df <- exon.resid.df[, isMapped.vector]
# merge the two dataframe by matching gene.df$EGO and rownames(exon.df)
merged.resid.df <- exon.resid.subset.df %>%
  mutate(EGO = as.integer(rownames(.))) %>%
  inner_join(gene.resid.df, by = "EGO") %>%
  column_to_rownames(var = "EGO")

# Spread genes into 100 output files and put relevant bins in the same output file as the gene.
total_genes <- length(allGenesCoveredByanyBins)
genes_per_file <- ceiling(total_genes / numSlice)

# Create a list to store gene groups
gene_groups <- split(allGenesCoveredByanyBins, 
                     rep(1:ceiling(total_genes/genes_per_file), each = genes_per_file, length.out = total_genes))
# For each group of genes length(gene_groups)
for (i in 1:length(gene_groups)) {
  
  # List to hold dataframes for each gene group
  dfs <- list()
  
  for (gene in gene_groups[[i]]) {
    # Identify columns that belong to the current gene or its exons
    matching_cols <- grepl(paste0("^", gene), colnames(merged.resid.df))
    
    # Slice the dataframe to include only the matching columns
    merged.resid.per.sigGene.df <- merged.resid.df[, matching_cols]
    
    # Add dataframe to list
    dfs[[gene]] <- merged.resid.per.sigGene.df
  }
  
  # Save the combined slice for conditional Genesis
  saveRDS(dfs, file = file.path(output.conditional.genesis.dir, paste0("conditionalSlice_", sprintf("%03d", i), ".RDS")))
}
write(paste("----------------END ", trait, "-------------------"), logPath, append=TRUE)