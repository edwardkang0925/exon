library(tidyverse)
  
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

# TODO: take this as command line argument and write slurm script?
trait = "lnTG"
FDR = 0.01
logPath <- "outputs/compareGeneBin.txt"


traits <- c("lnTG", "adjTC", "BMI", "fhshdl", "fvc") # traits with sig genes
for(trait in traits){
  write(paste0("----------",trait,"------------"), logPath, append = TRUE)
  # output path
  output.conditional.genesis.dir <- file.path("outputs", "conditionalGenesisInput", trait)
  
  # necessary files (TWAS results)
  exon_level_summary_file_path <- file.path("outputs", "combined", trait, paste0("combinedGenesis_", trait, ".csv"))
  gene_level_summary_file_path <- file.path("data", "geneLvMMAPresults", 
                                    paste0("parsed_output_twas_input_visit_1_adjusted_", trait,
                                           "_z_scores_p_values_corrected_significant_genes_bonferroni.csv"))
  exon_level_residual_file_path <- file.path("outputs", "afterRegression", trait, "combined", "residualExonExpression.RDS")
  gene_level_residual_file_path <- file.path("data", "MMAP_INPUT_FINAL_ROUND", paste0("twas_input_visit_1_adjusted_",trait,".csv"))
  # Read summary data of a trait
  exon.pvals.df <- read.csv(exon_level_summary_file_path)
  gene.pvals.df <- read.csv(gene_level_summary_file_path)
  exon.resid.df <- readRDS(exon_level_residual_file_path)
  gene.resid.df <- read.csv(gene_level_residual_file_path)
  
  # bin-level significance threshold
  exon.level.sig.threshold <- FDR / dim(exon.resid.df)[2]
  
  # collect significant genes from gene-level TWAS 
  # TODO: ENSG00000185291.12 is duplicated for lnTG. probably, PAR_Y thing. check this in exon
  # For now, using unique(), print number of sig bins 
  sigGenes <- unique(gene.pvals.df$Genes)
  sigExon.df <- exon.pvals.df[exon.pvals.df$pval_bacon < exon.level.sig.threshold,]
  write(paste(trait, "has", length(sigGenes), "significant genes"), logPath, append = TRUE)
  write(paste(trait, "has", dim(sigExon.df)[1], "significant bins"), logPath, append = TRUE)
  
  
  # print GIF 
  GIF_nominal <- inflation(exon.pvals.df$pval)
  write(paste(trait, "bin-level GIF (nominal) =", GIF_nominal), logPath, append = TRUE)
  
  GIF_bacon <- inflation(exon.pvals.df$pval_bacon)
  write(paste(trait, "bin-level GIF (bacon) =", GIF_bacon), logPath, append = TRUE)
  
  # Check exon -> gene coverage & sig-gene coverage
  sigExon.df <- sigExon.df <- sigExon.df %>%
    mutate(gene = mapBinToGene(exon))
  genes.containing.sigExons <- unique(sigExon.df$gene)
  write(paste(trait, "#unique genes containing sig bin:", length(genes.containing.sigExons)), logPath, append=TRUE)
  genesCoveredByBins <- intersect(genes.containing.sigExons, sigGenes) # we will output this number of slice for conditional analysis
  write(paste(trait, "#unique significant genes containing sig bin:", length(genesCoveredByBins)), logPath, append=TRUE)
  
  ## merge summary df and sort by name to group gene-bins.
  #first subset relevent slice from both residual dfs including the subject_ID (EGO)
  gene.resid.subset.df <- gene.resid.df[, union("EGO", colnames(gene.resid.df)[colnames(gene.resid.df) %in% genesCoveredByBins])]
  # Subsetting exon residual file: Apply the function to the column names of exon.resid.df
  mapped_gene_ids <- sapply(colnames(exon.resid.df), mapBinToGene)
  # Find which of these gene IDs are in genesCoveredByBins
  cols_to_keep <- mapped_gene_ids %in% genesCoveredByBins
  if(sum(cols_to_keep) == 0){
    write(paste(trait, "TERMINATE preparing conditional genesis as there is no significant gene"), logPath, append=TRUE)
    next
  }
  # create output dirs if not existing. this block is placed in the middle of the script instead of the top because if there is no sigBin, no reason to create directory.
  if(!dir.exists(output.conditional.genesis.dir)){
    dir.create(output.conditional.genesis.dir, recursive=TRUE)
  }
  
  # Subset exon.resid.df to only include these columns
  exon.resid.subset.df <- exon.resid.df[, cols_to_keep]
  # merge the two dataframe by matching gene.df$EGO and rownames(exon.df)
  merged.resid.df <- exon.resid.subset.df %>%
    mutate(EGO = as.integer(rownames(.))) %>%
    inner_join(gene.resid.subset.df, by = "EGO") %>%
    column_to_rownames(var = "EGO")
  
  ## Now split the merged dataframe per significant genes: rows = subject, columns = 1 sig gene and related bins.
  
  # For each significant genes containing significant bin
  for (gene in genesCoveredByBins) {
    
    # Identify columns that belong to the current gene or its exons
    matching_cols <- grepl(paste0("^", gene), colnames(merged.resid.df))
    
    # Slice the dataframe to include only the matching columns
    merged.resid.per.sigGene.df <- merged.resid.df[, matching_cols]
    
    # print how many bins exist within a sigGene containing sigBin.
    numBin <- dim(merged.resid.per.sigGene.df)[2] - 1 # subtract 1 because one is a gene not bin.
    numSigBin <- dim(sigExon.df[sigExon.df$gene == gene,])[1]
    if(numSigBin <= 3){ # print sigbinID if there are only couple of them.
      sigBins <- sigExon.df[sigExon.df$gene == gene, "exon"] # CHECK if this works
      write(paste(trait, "sigGeneContainingSigBin:", gene, "numBins:", numBin,
                  "numSigBins:", numSigBin, paste(sigBins, collapse = " ")), logPath, append=TRUE)
    } else{
      write(paste(trait, "sigGeneContainingSigBin:", gene, "numBins:", numBin,
                  "numSigBins:", numSigBin), logPath, append=TRUE)
    }
    
    # Save the slice for conditional Genesis
    saveRDS(merged.resid.per.sigGene.df, file = file.path(output.conditional.genesis.dir, paste0(gene, ".RDS")))
  }
  write(paste("----------------END ", trait, "-------------------"), logPath, append=TRUE)
}  
