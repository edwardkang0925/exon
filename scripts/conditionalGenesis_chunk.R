.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
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
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")
parser <- add_option(parser, c("--exon.path"), type="character",
                     default="outputs/conditionalGenesisInput_all/BMI/conditionalSlice_001.RDS",
                     help="path to exon slice",
                     metavar="'/path/to/file/conditionalSlice_###.RDS'")
parser <- add_option(parser, c("--outputFile"), type="character",
                     default="outputs/conditionalGenesisOut_all/BMI/results_001.RDS",
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
output.path <- opt$outputFile
trait <- opt$trait # New argument for the trait

# create output dir
directory_path <- dirname(output.path)
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

# path to required files: 
kinship.matrix.path <- paste0("/project/data/kinship/kinship_",trait,".RDS")
kinship.matrix <- readRDS(kinship.matrix.path)
trait.df.path <- paste0("/project/data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_", trait,".csv")
trait.df <- read.csv(trait.df.path)
print(kinship.matrix)


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
  
  # Loop through each exon: do trait ~ exon | gene considering kinship_matrix
  for (exon_col in exon_columns) {
    model_data <- data.frame(
      trait = merged_df[[trait]],
      exon = merged_df[[exon_col]],
      gene = merged_df[[gene_column]],
      row.names = rownames(merged_df)
    )
    tryCatch({
      # Fit the null model
      null_model <- fitNullModel(
        x = model_data,
        outcome = "trait",
        covars = "gene",
        cov.mat = kinship.matrix
      )
      
      # Fit the alternative model and perform the association test
      alt_model <- fitNullModel(
        x = model_data,
        outcome = "trait",
        covars = c("gene", "exon"),
        cov.mat = kinship.matrix
      )
      
      # Calculate LRT
      lrt <- 2 * (alt_model$logLik - null_model$logLik)
      p_value_lrt <- pchisq(lrt, df = 1, lower.tail = FALSE)
      p_value_lrt_bon <- min(1.0, p_value_lrt * length(exon_columns))
      
      # Append the results to the results data frame
      results <- rbind(results, data.frame(exon = exon_col, p_lrt = p_value_lrt,
                                          p_lrt_bon = p_value_lrt_bon, lrt_stat=lrt))
    }, error = function(e) {
      print(paste0("Error: ", e))
      print(paste0("Exon: ", exon_col))
    })
  }
  
  # Add the results to the results list
  results_list[[i]] <- results
}

# Combine all results across genes within the input slice into a single data frame 
combined_results <- do.call(rbind, results_list)
print(dim(combined_results))
# Save the combined results
saveRDS(combined_results, file = output.path)
