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
                     default="outputs/geneBinCorrelation/slices/results_001.RDS",
                     help="output filename including path",
                     metavar="/path/to/output/outputfilename")


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

# create output dir
directory_path <- dirname(output.path)
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

# Load the exon expression data
exon_data <- readRDS(exon.path)

# Initialize an empty results list
results_list <- list()

# Iterate over each element (data frame) in the exon_data list
for (i in seq_along(exon_data)) {
  # Get the current data frame from the list
  current_df <- exon_data[[i]]
  
  # Find columns to work with
  colnames_data <- colnames(current_df)
  gene_column <- grep("^ENSG\\d+\\.\\d+$", colnames_data, value = TRUE)
  exon_columns <- grep("^ENSG\\d+\\.\\d+\\.\\d+$", colnames_data, value = TRUE)
  
  # Initialize an empty results data frame
  results <- data.frame()
  
  # Loop through each exon: do trait ~ exon | gene considering kinship_matrix
  for (exon_col in exon_columns) {
    # calculate correlation between gene and bin. 
    correlation <- cor(current_df[gene_column], current_df[exon_col], method = "pearson")
    
    
    # Append the results to the results data frame
    results <- rbind(results, data.frame(exon = exon_col, cor = correlation[[1]]))
  }
  
  # Add the results to the results list
  results_list[[i]] <- results
}

# Combine all results across genes within the input slice into a single data frame 
combined_results <- do.call(rbind, results_list)
print(dim(combined_results))
# Save the combined results
saveRDS(combined_results, file = output.path)
