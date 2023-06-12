.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript
library(dplyr)
library(tibble)
library(stats)
library(optparse)

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
parser <- OptionParser()
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")
parser <- add_option(parser, c("--numSlice"), type="numeric", default=100,
                      help="number of slices")
parser <- add_option(parser, c("--pvalThreshold"), type="numeric", default=0.01,
                      help="pre-bonferroni pvalThreshold for bin-level significance")
parser <- add_option(parser, c("--mergedResidualsPath"), type="character", default="",
                      help="path to merged residuals")
parser <- add_option(parser, c("--outputDir"), type="character", default="outputs/conditionalGenesisInput/BMI/",
                      help="output directory")
opt = parse_args(parser)

# Assign the first argument to trait
trait <- opt$trait
numSlice = opt$numSlice
FDR = opt$pvalThreshold

# output path
output.conditional.genesis.dir <- opt$outputDir

# create output dirs if not existing. 
if(!dir.exists(output.conditional.genesis.dir)){
  dir.create(output.conditional.genesis.dir, recursive=TRUE)
}

# necessary files (TWAS results)
exon_level_residual_file_path <- opt$mergedResidualsPath
gene_level_residual_file_path <- file.path("/project/data", "MMAP_INPUT_FINAL_ROUND", paste0("twas_input_visit_1_adjusted_",trait,".csv"))
# Read summary data of a trait SLOW
exon.resid.df <- readRDS(exon_level_residual_file_path)
gene.resid.df <- read.csv(gene_level_residual_file_path)
allGenes <- unique(colnames(gene.resid.df)[grepl("ENSG", colnames(gene.resid.df))])

# identify all genes containing any bin. (preparing input for conditional analysis on all bins within a gene)
genesWithBin <- unique(lapply(colnames(exon.resid.df), mapBinToGene))
allGenesCoveredByanyBins <- intersect(allGenes, genesWithBin)

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
