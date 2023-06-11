# Obtain the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign the first argument to trait
trait <- args[1]

# split_exon_files.R
exon.path <- paste0("outputs/afterRegression/", trait, "/combined/residualExonExpression.RDS")
phenotype.df.path <- paste0("data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_", trait,".csv")
dir_path <- paste0("outputs/genesisInput/",trait,"/")

# Load the exon expression data and phenotype df
exon_data <- readRDS(exon.path)
phenotype.df <- read.csv(phenotype.df.path)
# Identify the column names that are exon ids
exon_cols <- grep("ENSG", colnames(exon_data), value = TRUE)
print(paste0("Dim of combined residual expression: ", dim(exon_data)))
# Calculate the number of exons per file
exons_per_file <- ceiling(length(exon_cols) / 100)

# Split the exons into 200 groups
split_exons <- split(exon_cols, ceiling(seq_along(exon_cols) / exons_per_file))

# Check if the directory exists
if (!dir.exists(dir_path)) {
  # If it does not exist, create it
  dir.create(dir_path)
}


# Create each file appending trait value
for(i in seq_along(split_exons)) {
  file_data <- exon_data[, split_exons[[i]]]
  # Convert rownames to a column
  file_data$subject <- rownames(file_data)
  # Merge the two dataframes
  merged_df <- merge(file_data, phenotype.df, by = "subject")
  # Assign rownames before dropping the column
  rownames(merged_df) <- merged_df$subject
  # drop the column `subject`
  merged_df$subject <- NULL
  saveRDS(merged_df, sprintf("outputs/genesisInput/%s/exon_data_%03d.RDS", trait, i))
}

