library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)  # For the gather function used later
library(GENESIS)

# Investigating bins with low pvalue.
# for such bins, plot bin expression vs gene expression plot to see what is going on.

# CHANGE correlation threshold (and inequality sign) and the output directory accordingly!!

outputPath = "tmpdata/exonWideSigBins/"

## Get list of bins with high correlation and low conditional pvalue.
df_corr = readRDS("data/combinedCorrelation_BMI.RDS")
df_corr_high = df_corr #[df_corr$cor > 0.95] # 403 bins for > 0.95 CHANGE ME

df_lrt_p = readRDS("data/combinedConditionalGenesis_BMI.RDS")

# filter bins with p_lrt_bon < 0.05
df_lrt_lowp = df_lrt_p[df_lrt_p$p_lrt < 0.05 / 586385,]

# find bins with high corr and low p
bins_highCor_and_lowP <- df_lrt_lowp$exon

# read bin expression data
binExpressionFilePath = "data/residualExonExpression.RDS"
df_binExpression = readRDS(binExpressionFilePath)

# read trait data
traitPath = "data/adjusted_BMI.csv"
df_trait = read.csv(traitPath)

# Assuming your df_trait dataframe has two columns, subject and BMI,
# and the row names of df_binExpression match with the subject column of df_trait.
bmi_values <- df_trait$BMI[match(rownames(df_binExpression), df_trait$subject)]

# read gene expression data
trait <- "BMI" # I added this since it seems 'trait' was not defined in the provided code
gene_level_residual_file_path <- file.path("data", paste0("twas_input_visit_1_adjusted_", trait, ".csv"))
df_geneExpression = read.csv(gene_level_residual_file_path)

# path to required files:
kinship.matrix.path <- paste0("data/kinship_",trait,".RDS")
kinship_matrix <- readRDS(kinship.matrix.path)
kinship_ids <- rownames(kinship_matrix)
names(kinship_ids) <- kinship_ids

# Loop through each bin of interest
for(bin in bins_highCor_and_lowP) {

  # Subset bin expression for the current bin
  bin_expression <- df_binExpression[[bin]]

  # Convert bin name to gene name by removing the last ".X" part
  gene_name <- unlist(strsplit(bin, "\\."))[1:(length(unlist(strsplit(bin, "\\."))) - 1)] %>% paste(collapse = ".")

  # Extract gene expression for the current gene
  gene_expression <- df_geneExpression[[gene_name]]

  # Create dataframe for the current bin-gene pair
  df_current <- data.frame(exonExpression = bin_expression, geneExpression = gene_expression, BMI = bmi_values, id = df_geneExpression$EGO)
  rownames(df_current) <- df_current$id

  # 1. Scatter plot of gene expression vs. trait with the bin expression overlayed
  p1 <- ggplot(df_current, aes(x = geneExpression, y = BMI, color = exonExpression)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") + # Regression line added here
    labs(title = paste0("Gene vs. Trait colored by ", bin), x = "Gene Expression", y = "BMI") +
    scale_color_gradient(name = "Bin Expression", low = "blue", high = "red")
  ggsave(filename = file.path(outputPath, paste0("1/gene_vs_trait_colored_by_", bin, ".png")), plot = p1, width = 6, height = 4)

  # 2. Scatter plot of the residual bin expression vs. the trait
  # First, regress out the gene expression from the bin expression

  # ensure the id of kinship matrix and expression dataframe matches.
  kinship_row_indices <- kinship_ids[as.character(df_current$id)]
  kinship_matrix_ordered <- kinship_matrix[kinship_row_indices, kinship_row_indices]

  # regress out gene expression on bin expression
  null_model <- fitNullModel(
    x = df_current,
    outcome = "exonExpression",
    covars = "geneExpression",
    cov.mat = kinship_matrix_ordered
  )

  residual_bin_expression <- null_model$fit$resid.marginal
  df_current$residual_bin_expression = residual_bin_expression

  p2 <- ggplot(df_current, aes(x = residual_bin_expression, y = BMI)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") + # Regression line added here
    theme_minimal() +
    labs(title = paste0("Residual Bin Expression vs. Trait for ", bin), x = "Residual Bin Expression", y = "BMI")
  ggsave(filename = file.path(outputPath, paste0("2/residual_bin_expression_vs_trait_", bin, ".png")), plot = p2, width = 6, height = 4)

  # 3. Scatter plot of residual bin expression vs. residual BMI
  # First, regress out the bin expression from the BMI to get residuals
  null_model2 <- fitNullModel(
    x = df_current,
    outcome = "BMI",
    covars = "geneExpression",
    cov.mat = kinship_matrix_ordered
  )

  residual_bmi = null_model2$fit$resid.marginal
  df_current$residual_BMI = residual_bmi # Adding it to df_current for plotting

  p3 <- ggplot(df_current, aes(x = residual_bin_expression, y = residual_BMI)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") + # Regression line added here
    theme_minimal() +
    labs(title = paste0("Residual Bin Expression vs. Residual BMI for ", bin), x = "Residual Bin Expression", y = "Residual BMI")
  ggsave(filename = file.path(outputPath, paste0("3/residual_bin_expression_vs_residual_BMI_", bin, ".png")), plot = p3, width = 6, height = 4)
}


