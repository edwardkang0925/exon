# Load necessary library
library(data.table)
library(ggplot2)
library(dplyr)

# Investigating bins with low pvalue. 
# for such bins, plot bin expression vs gene expression plot to see what is going on.

## Get list of bins with high correlation and low conditional pvalue. 
df_corr = readRDS("outputs/geneBinCorrelation/combinedCorrelation_BMI.RDS")
df_corr_high = df_corr[df_corr$cor > 0.95,] # 403 bins
head(df_corr_high) # 2 columns: [exon, cor]
df_lrt_p = readRDS("outputs/combinedConditionalGenesis/BMI/combinedConditionalGenesis_BMI.RDS")
head(df_lrt_p) # 4 columns: [exon, p_lrt, p_lrt_bon, lrt_stat]

# filter bins with p_lrt_bon < 0.05, currently 3471 bins
df_lrt_lowp = df_lrt_p[df_lrt_p$p_lrt_bon < 0.05,] 
dim(df_lrt_lowp)

# find bins with high corr and low p
bins_highCor_and_lowP <- intersect(df_corr_high$exon , df_lrt_lowp$exon)
length(bins_highCor_and_lowP) 

# 26 bins with high corr and low p. 
## For each bin, grab a vector of bin expression and gene expression across individuals

# read bin expression data
binExpressionFilePath = "outputs/afterRegression/BMI/combined/residualExonExpression.RDS"
df_binExpression = readRDS(binExpressionFilePath)

# read gene expression data
gene_level_residual_file_path <- file.path("data", "MMAP_INPUT_FINAL_ROUND", paste0("twas_input_visit_1_adjusted_",trait,".csv"))
df_geneExpression = read.csv(gene_level_residual_file_path)

# grab vectors, output 26 files each having two columns [exonExpression, geneExpression]
# Loop through each bin of interest
for(bin in bins_highCor_and_lowP) {
  
  # Subset bin expression for the current bin
  bin_expression <- df_binExpression[[bin]]
  
  # Convert bin name to gene name by removing the last ".X" part
  gene_name <- unlist(strsplit(bin, "\\."))[1:(length(unlist(strsplit(bin, "\\."))) - 1)] %>% paste(collapse = ".")
  
  # Extract gene expression for the current gene
  gene_expression <- df_geneExpression[[gene_name]]
  
  # Create dataframe for the current bin-gene pair
  df_current <- data.frame(exonExpression = bin_expression, geneExpression = gene_expression)
  
  # Save to a CSV file
  write.csv(df_current, file.path("tmpdata/highCorrLowPbins/", paste0("expression_data_", bin, ".csv")), row.names = FALSE)
  
  # Create scatter plot
  plot <- ggplot(df_current, aes(x = geneExpression, y = exonExpression)) + 
    geom_point(alpha = 0.5) +
    theme_minimal() + 
    labs(title = paste0("Scatter plot for ", bin), x = "Gene Expression", y = "Bin Expression")
  
  # Save the plot to a file
  ggsave(filename = file.path("tmpdata/highCorrLowPbinsPlots/", paste0("scatter_plot_", bin, ".png")), plot = plot, width = 6, height = 4)
}
