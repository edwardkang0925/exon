#!/usr/bin/env Rscript
library(optparse)
library(SeqArray)
library(data.table)
#library(kinship2)
library(dplyr)
library(stringr)
library(tibble)
library(foreach)
library(GENESIS)

parser <- OptionParser()
parser <- add_option(parser, c("--exon.path"), type="character",
                     default="outputs/conditionalGenesisInput_all/BMI/conditionalSlice_001.RDS",
                     help="path to exon slice",
                     metavar="'/path/to/file/conditionalSlice_###.RDS'")
parser <- add_option(parser, c("--output.path"), type="character",
                     default="outputs/conditionalGenesisOut_regression_all/BMI/results_001.RDS",
                     help="output filename including path",
                     metavar="/path/to/output/outputfilename")
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")

opt = parse_args(parser)

# check cmd line arguments
x = foreach(
  i = names(opt)
) %do% {
  input_value=opt[[i]]
  if(input_value==''){
    stop(sprintf("ARGUMENT --%s IS REQUIRED",i))
  }else if(i %in% c('dna','rna')){
    if(!file.exists(input_value)){
      stop(sprintf("FILE %s DOES NOT EXIST",input_value))
    }
  }
}

# read command line input
exon.path <- opt$exon.path 
output.path <- opt$output.path
trait <- opt$trait # New argument for the trait

# create output dir
directory_path <- dirname(output.path)
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

# path to required files: 
kinship.matrix.path <- paste0("data/kinship/kinship_",trait,".RDS")
kinship_matrix <- readRDS(kinship.matrix.path)
trait.df.path <- paste0("data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_", trait,".csv")
trait.df <- read.csv(trait.df.path)


# Load the exon expression data
exon_data <- readRDS(exon.path)

# Initialize an empty results list
results_list <- list()

# Iterate over each element (data frame) in the exon_data list
for (i in seq_along(exon_data)) {
  # Get the current data frame from the list
  current_df <- exon_data[[i]]
  
  # Merge with the trait data frame
  current_df$subject <- rownames(current_df)
  merged_df <- merge(current_df, trait.df, by="subject") %>%
    column_to_rownames(var="subject")
  
  
  # Find columns to work with
  colnames_data <- colnames(merged_df)
  gene_column <- grep("^ENSG\\d+\\.\\d+$", colnames_data, value = TRUE)
  exon_columns <- grep("^ENSG\\d+\\.\\d+\\.\\d+$", colnames_data, value = TRUE)
  
  # Initialize an empty results data frame
  results <- data.frame()
  
  # Initialize data for the first regression (trait ~ gene)
  model_data <- data.frame(
    trait = merged_df[[trait]],
    gene = merged_df[[gene_column]],
    row.names = rownames(merged_df)
  )
  
  # First regression: regress out gene's effect and calculate the conditional residuals
  regression_gene_model <- fitNullModel(
    x = model_data,
    outcome = "trait",
    covars = "gene",
    cov.mat = kinship_matrix
  )
  gene_residuals <- regression_gene_model$fit['resid.conditional']
  
  # Loop through each exon: do trait ~ exon | gene considering kinship_matrix
  for (exon_col in exon_columns) {
    model_data <- data.frame(
      resid.conditional = gene_residuals,
      exon = merged_df[[exon_col]],
      row.names = rownames(merged_df)
    )
    
    # Second regression: residual (after regressing out gene)  ~ bin 
    regression_bin_model <- fitNullModel(
      x = model_data,
      outcome = "resid.conditional",
      covars = "exon",
      cov.mat = kinship_matrix
    )
    beta <- regression_bin_model$fixef["exon", "Est"]
    se <- regression_bin_model$fixef["exon", "SE"]
    pval <- regression_bin_model$fixef["exon", "pval"]
    
    # Append the results to the results data frame
    results <- rbind(results, data.frame(exon = exon_col, p_seq_regression = pval,
                                         beta = beta, se = se))
  }
  
  # Add the results to the results list
  results_list[[i]] <- results
}

# Combine all results across genes within the input slice into a single data frame 
combined_results <- do.call(rbind, results_list)
print(dim(combined_results))
# Save the combined results
saveRDS(combined_results, file = output.path)
