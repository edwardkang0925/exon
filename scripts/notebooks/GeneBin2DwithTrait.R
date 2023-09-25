# Load necessary libraries
library(data.table)
library(ggplot2)
library(dplyr)

# Investigating bins with low pvalue. 
# for such bins, plot bin expression vs gene expression plot to see what is going on.

# CHANGE correlation threshold (and inequality sign) and the output directory accordingly!!

outputPath = "tmpdata/medCorrLowPbinsPlots/"

## Get list of bins with high correlation and low conditional pvalue. 
df_corr = readRDS("outputs/geneBinCorrelation/combinedCorrelation_BMI.RDS")
df_corr_high = df_corr[df_corr$cor >= 0.1 & df_corr$cor <= 0.5, ] # 403 bins for > 0.95 CHANGE ME 

df_lrt_p = readRDS("outputs/combinedConditionalGenesis/BMI/combinedConditionalGenesis_BMI.RDS")

# filter bins with p_lrt_bon < 0.05
df_lrt_lowp = df_lrt_p[df_lrt_p$p_lrt_bon < 0.05,] 

# find bins with high corr and low p
bins_highCor_and_lowP <- intersect(df_corr_high$exon, df_lrt_lowp$exon)

# read bin expression data
binExpressionFilePath = "outputs/afterRegression/BMI/combined/residualExonExpression.RDS"
df_binExpression = readRDS(binExpressionFilePath)

# read trait data
traitPath = "data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_BMI.csv"
df_trait = read.csv(traitPath)

# Assuming your df_trait dataframe has two columns, subject and BMI, 
# and the row names of df_binExpression match with the subject column of df_trait.
bmi_values <- df_trait$BMI[match(rownames(df_binExpression), df_trait$subject)]

# read gene expression data
trait <- "BMI" # I added this since it seems 'trait' was not defined in the provided code
gene_level_residual_file_path <- file.path("data", "MMAP_INPUT_FINAL_ROUND", paste0("twas_input_visit_1_adjusted_", trait, ".csv"))
df_geneExpression = read.csv(gene_level_residual_file_path)

# Loop through each bin of interest
for(bin in bins_highCor_and_lowP) {
  
  # Subset bin expression for the current bin
  bin_expression <- df_binExpression[[bin]]
  
  # Convert bin name to gene name by removing the last ".X" part
  gene_name <- unlist(strsplit(bin, "\\."))[1:(length(unlist(strsplit(bin, "\\."))) - 1)] %>% paste(collapse = ".")
  
  # Extract gene expression for the current gene
  gene_expression <- df_geneExpression[[gene_name]]
  
  # Create dataframe for the current bin-gene pair
  df_current <- data.frame(exonExpression = bin_expression, geneExpression = gene_expression, BMI = bmi_values)
  
  # Create 2D scatter plot with color representing BMI
  plot <- ggplot(df_current, aes(x = geneExpression, y = exonExpression, color = BMI)) + 
    geom_point(alpha = 0.5) +
    theme_minimal() + 
    labs(title = paste0("Scatter plot for ", bin), x = "Gene Expression", y = "Bin Expression") +
    scale_color_gradient(name = "BMI", low = "blue", high = "red")  # Adjust the colors as needed
  
  # Save the plot to a file
  ggsave(filename = file.path(outputPath, paste0("color_scatter_plot_", bin, ".png")), plot = plot, width = 6, height = 4)
}
