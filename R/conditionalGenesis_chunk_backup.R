.libPaths(c("/home/edwardkang/R/x86_64-conda-linux-gnu-library/4.2", .libPaths()))
library(GENESIS)
library(SeqArray)
library(data.table)
#library(kinship2)
library(dplyr)
library(stringr)
library(tibble)

# genesis_analysis.R
args <- commandArgs(trailingOnly = TRUE)
exon.path <- args[1]
output.path <- args[2]
trait <- args[3] # New argument for the trait

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

# Initialize an empty results data frame
results <- data.frame()

# check if trait is in the file: FIXME: when preparing input, I didn't attach trait value, so appending in the forloop below.
#print(trait %in% colnames(exon_data))
exon_data$subject <- rownames(exon_data)
exon_data <- merge(exon_data, trait.df, by="subject") %>%
  column_to_rownames(var="subject")

# Find columns to work with
colnames_data <- colnames(exon_data)
gene_column <- grep("^ENSG\\d+\\.\\d+$", colnames_data, value = TRUE)
exon_columns <- grep("^ENSG\\d+\\.\\d+\\.\\d+$",colnames_data, value = TRUE)

# Loop through each exon: do trait ~ exon | gene considering kinship_matrix
for (exon_col in exon_columns) {
  model_data <- data.frame(trait = exon_data[[trait]],  
                           exon = exon_data[[exon_col]],
                           gene = exon_data[[gene_column]],
                           row.names = rownames(exon_data))
  # Fit the null model
  null_model <- fitNullModel(
    x = model_data,
    outcome = "trait",
    covars = "gene",  
    cov.mat = kinship_matrix
  )
  
  # Fit the alternative model and perform the association test
  alt_model <- fitNullModel(
    x = model_data,
    outcome = "trait",
    covars = c("gene","exon"), 
    cov.mat = kinship_matrix
  )
  
  # calculate LRT
  lrt <- 2*(alt_model$logLik- null_model$logLik)
  p_value_lrt <- pchisq(2 * lrt, df = 1, lower.tail = FALSE)
  p_value_lrt_bon <- min(1.0, p_value_lrt * length(exon_columns))
  
  # Append the results to the results data frame
  results <- rbind(results, data.frame(exon = exon_col, p_lrt = p_value_lrt, p_lrt_bon = p_value_lrt_bon))
}


# Save the results
saveRDS(results, file = output.path)
