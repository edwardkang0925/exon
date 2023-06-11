library(GENESIS)
library(SeqArray)
library(data.table)
library(kinship2)


# path to required files.
exon.path <- "expected/input/twas_input_visit_1_adjusted_fhshdl.csv"
pedigree.path <- "data/mmap_input/mmap.ped.v5.csv"

# create a kinship matrix from the pedigree file.
pedigree_data <- read.csv(pedigree.path)
ped <- with(pedigree_data, pedigree(id = ID, dadid = FA, momid = MO, sex = SEX, famid = FAMID))
kinship_matrix <- kinship(ped)

# Load the exon expression data
exon_data <- fread(exon.path)
not_ENSG_cols <- grep("ENSG", colnames(exon_data), invert = TRUE, value = TRUE)
not_ENSG_cols

# Set rownames as subject ID
rownames(exon_data) <- exon_data$EGO

# Remove the individual column
if("EGO" %in% colnames(exon_data)) {
  exon_data$EGO <- NULL
}

# subset kinship matrix to match the exon_data
# Get the common ids
common_ids <- intersect(rownames(exon_data), rownames(kinship_matrix))
# Subset the kinship matrix
kinship_matrix <- kinship_matrix[common_ids, common_ids]
# match the order of rownames between kinshipmatrix and exon_data
kinship_matrix <- kinship_matrix[rownames(exon_data), rownames(exon_data)]


# Initialize an empty results data frame
results <- data.frame()

# Loop through each exon and do association test: NOTE: the p-val is not corrected for the multiple testing. 
for (exon in colnames(exon_data)) {
  print(exon)
  if(grepl("ENSG", exon)){
    # Prepare data for the model
    model_data <- data.frame(fhshdl = exon_data[["fhshdl"]], 
                             exon = exon_data[[exon]], 
                             row.names = rownames(exon_data))
    
    # Fit the null model
    null_model <- fitNullModel(
      x = model_data,
      outcome = "fhshdl",
      covars = NULL,  # No covariates except the intercept
      cov.mat = kinship_matrix
    )
    
    # Fit the full model with exon as a covariate
    full_model <- fitNullModel(
      x = model_data,
      outcome = "fhshdl",
      covars = "exon",
      cov.mat = kinship_matrix
    )
    # Perform a likelihood ratio test
    lrt <- 2*(full_model$logLik- null_model$logLik)
    print(lrt)  # Print out the lrt value
    p_value <- pchisq(2 * lrt, df = 1, lower.tail = FALSE)
    # Append the results to the results data frame
    results <- rbind(results, data.frame(exon = exon, p_value = p_value))
  }
}
# View the results
results
