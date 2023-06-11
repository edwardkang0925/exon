# refactoring genesis_chunk.R 05.22.2023 to use preprocessed kinship matrix. Done
library(GENESIS)
library(SeqArray)
library(data.table)
library(kinship2)

# genesis_analysis.R
args <- commandArgs(trailingOnly = TRUE)
exon.path <- args[1]
output.path <- args[2]
trait <- args[3] # New argument for the trait

# delete later
trait <- "fvc"
exon.path <- paste0("outputs/genesisInput/",trait,"/exon_data_001.RDS")
# output.path <- paste0("outputs/genesisOut/",trait,"/results_001.RDS")


# path to required files: TODO: for each trait, preprocess the kinship matrix to reduce the computation burden.
kinship.matrix.full.path <- paste0("data/kinship/kinship_",trait,".RDS")
kinship_matrix <- readRDS(kinship.matrix.full.path)

# Load the exon expression data
exon_data <- readRDS(exon.path)
not_ENSG_cols <- grep("ENSG", colnames(exon_data), invert = TRUE, value = TRUE)
not_ENSG_cols

# Initialize an empty results data frame
results <- data.frame()

# check if trait is in the file
print(trait %in% colnames(exon_data)) 

# Loop through each exon and do association test
for (exon in colnames(exon_data)) {
  if(grepl("ENSG", exon)){
    model_data <- data.frame(trait = exon_data[[trait]],  
                             exon = exon_data[[exon]], 
                             row.names = rownames(exon_data))
    print(sum(is.na(model_data)))
    
    full_model <- fitNullModel(
      x = model_data,
      outcome = "trait",  # This looks for a column in `x` that is named with the value of `trait` (ex. ABI)
      covars = "exon", # This looks for a column in `x` named "exon"
      cov.mat = kinship_matrix
    )
    
    beta <- full_model$fixef["exon", "Est"]
    se <- full_model$fixef["exon", "SE"]
    pval <- full_model$fixef["exon", "pval"]
    
    results <- rbind(results, data.frame(exon = exon, beta = beta, se = se, pval = pval))
  }
}

# # create output dir
# if (!dir.exists(output.path)) {
#   dir.create(output.path, recursive = TRUE)
# }
# 
# # Save the results
# saveRDS(results, file = output.path)
