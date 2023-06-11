# Load necessary library
library(data.table)
library(bacon)
library(tidyverse)
library(ggplot2)


# it is validating whether Genesis and MMAP produce comparable results. 

# Specify the path to the directory containing the RDS files
path_to_rds_files <- "expected/genesisOut/"
outputPATH <- "expected/combined/genesisCombinedGeneLevel.csv"
trait = "fhshdl"
hist.logbaconp.path <- paste0("outputs/plots/hist_genesis_mmap_logp_diff_",trait,".png")
hist.logbaconp.zoomed.path <- paste0("outputs/plots/hist_genesis_mmap_logp_diff_",trait,"_zoomed.png")


# Get a list of all RDS files in the directory
rds_files <- list.files(path = path_to_rds_files, pattern = "*.RDS", full.names = TRUE)

# Initialize an empty data table
all_data <- data.table()

# Loop over each RDS file
for (file in rds_files) {
  # Load the data from the RDS file
  data <- readRDS(file)
  
  # rbind the data to the existing data table
  all_data <- rbindlist(list(all_data, data))
}
# Calculate the z-score
#all_data$z_score <- all_data$beta / all_data$se

# Calculate the p-value from the z-score
#all_data$p_value_bs <- 2 * (1 - pnorm(abs(all_data$z_score)))


# save CSV file
write.csv(all_data, file = outputPATH)


# compare with the Sandeep generated file 
sandeep.path <- "expected/expectedOutput/parsed_output_twas_input_visit_1_adjusted_fhshdl.csv"
df.sandeep <- read.csv(sandeep.path)
# check for duplicated Genes name.
df.sandeep[duplicated(df.sandeep$Genes) | duplicated(df.sandeep$Genes, fromLast = TRUE), ]

# remove duplicated rows
df.sandeep <- df.sandeep[!duplicated(df.sandeep), ]
df.sandeep[duplicated(df.sandeep$Genes) | duplicated(df.sandeep$Genes, fromLast = TRUE), ]
# still have duplicated Genes because it has different p_vals. 

commonGenes <- intersect(all_data$exon, df.sandeep$Genes)
print(length(commonGenes))
commonGenes <- commonGenes[commonGenes != "ENSG00000182162.11"] # dropping duplicated gene for fair comparison
# since all_data_sorted has 14 more genes, subset to match the genes
all_data <- all_data[all_data$exon %in% commonGenes,]
df.sandeep <- df.sandeep[df.sandeep$Genes %in% commonGenes,]

print(dim(all_data))
print(dim(df.sandeep))

# sort by pvalue
all_data_sorted <- all_data %>% arrange(pval)
df_sandeep_sorted <- df.sandeep %>% arrange(p_vals)


# Count the number of matches
num_matches <- sum(all_data_sorted$exon == df_sandeep_sorted$Genes)
matching_index <- all_data_sorted$exon == df_sandeep_sorted$Genes
all_data_sorted[matching_index,]

# Print the number of matches
print(paste("Number of matches:", num_matches))


## compare logP vals before BACON
# Calculate -log10(p_value) for both dataframes
all_data_sorted$neg_log_pval_pre_bacon <- -log10(all_data_sorted$pval)
df_sandeep_sorted$neg_log_pval_pre_bacon <- -log10(df_sandeep_sorted$p_vals)

# Merge the two dataframes
names(df_sandeep_sorted)[names(df_sandeep_sorted) == "Genes"] <- "exon"
merged_df_pre_bacon <- merge(all_data_sorted, df_sandeep_sorted, by = "exon")

# Calculate the abs difference between -log10(p_value) of 'exon' and 'Genes'
merged_df_pre_bacon$abs_diff_neg_log_pval_pre_bacon <- abs(merged_df_pre_bacon$neg_log_pval_pre_bacon.x - merged_df_pre_bacon$neg_log_pval_pre_bacon.y)

# Calculate maximum difference
max_diff_pre_bacon <- max(merged_df_pre_bacon$abs_diff_neg_log_pval_pre_bacon, na.rm = TRUE)

# Plot the distribution of the difference
hist.plot.pre.bacon <- ggplot(merged_df_pre_bacon, aes(x = abs_diff_neg_log_pval_pre_bacon)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "white") +  # Adjust binwidth as needed
  geom_vline(aes(xintercept = max_diff_pre_bacon), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max_diff_pre_bacon, y = Inf, label = paste("Max Diff =", round(max_diff_pre_bacon, 2)), vjust = 2, hjust = 1.5, color = "red") +
  theme_minimal() +
  labs(x = "Absolute Difference in -log10(p-value)", y = "Count",
       title = paste0("GENESIS vs MMAP (", trait, "): abs diff -log10(P) pre-BACON"))
hist.plot.pre.bacon
ggsave(paste0("outputs/plots/hist_genesis_mmap_logp_diff_",trait,"_pre_bacon.png"), plot = hist.plot.pre.bacon, width = 10, height = 7, dpi = 300)

# plot (zoom-in version toward the extreme values)
cutoff = 0.2
filtered_df_pre_bacon <- merged_df_pre_bacon[merged_df_pre_bacon$abs_diff_neg_log_pval_pre_bacon >= cutoff,]

hist.zoom.plot.pre.bacon <- ggplot(filtered_df_pre_bacon, aes(x = abs_diff_neg_log_pval_pre_bacon)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "white") +  # Adjust binwidth as needed
  geom_vline(aes(xintercept = max_diff_pre_bacon), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max_diff_pre_bacon, y = Inf, label = paste("Max Diff =", round(max_diff_pre_bacon, 2)), vjust = 2, hjust = 1.5, color = "red") +
  theme_minimal() +
  labs(x = paste("Difference in -log10(P-value), cutoff=", cutoff), y = "Count",
       title = paste0("GENESIS vs MMAP (", trait, "): zoomed abs diff -log10(P) pre-bacon")) +
  coord_cartesian(xlim = c(0, max(filtered_df_pre_bacon$abs_diff_neg_log_pval)), ylim = c(0, 20))  # Set x limits and y limits

hist.zoom.plot.pre.bacon
ggsave(paste0("outputs/plots/hist_genesis_mmap_logp_diff_",trait,"_zoomed_pre_bacon.png"), plot = hist.zoom.plot.pre.bacon, width = 10, height = 7, dpi = 300)



colnames(merged_df_pre_bacon)

merged_df_pre_bacon_passing_mmap <- merged_df_pre_bacon[merged_df_pre_bacon$p_vals < 0.05/2*10^-4,]
merged_df_pre_bacon_passing_genesis <- merged_df_pre_bacon[merged_df_pre_bacon$pval < 0.05/2*10^-4,]

gene.intersect <- intersect(merged_df_pre_bacon_passing_mmap$exon, merged_df_pre_bacon_passing_genesis$exon)
length(gene.intersect)
merged_df_pre_bacon_passing_genesis[!merged_df_pre_bacon_passing_genesis$exon %in% gene.intersect,]



       
       
## BACON and compare logP vals.
# first create a bacon object:
mmap_z_score = qnorm(df_sandeep_sorted$p_vals/2, lower.tail=FALSE) *
  sign(df_sandeep_sorted$betas_genes / df_sandeep_sorted$se_genes)
genesis_z_score = qnorm(all_data_sorted$p_value/2, lower.tail=FALSE) *
  sign(all_data_sorted$beta / all_data_sorted$se)

mmap_bacon <- bacon(teststatistics = mmap_z_score, niter=10000L, nburnin = 3000L)
genesis_bacon <- bacon(teststatistics = genesis_z_score, niter=10000L, nburnin = 3000L)


# add BACON corrected pvalues
df_sandeep_sorted$adjusted_pvalues <- pval(mmap_bacon, corrected = TRUE)
all_data_sorted$adjusted_pvalues <- pval(genesis_bacon, corrected = TRUE)

# sort by adj-pvalue
all_data_sorted <- all_data_sorted %>% arrange(adjusted_pvalues)
df_sandeep_sorted <- df_sandeep_sorted %>% arrange(adjusted_pvalues)

# Count the number of matches
num_matches <- sum(all_data_sorted$exon == df_sandeep_sorted$Genes)
matching_index <- all_data_sorted$exon == df_sandeep_sorted$Genes

# Print the number of matches
print(paste("Number of matches:", num_matches))

# Print the total number of elements
total_elements <- length(all_data_sorted$exon)
print(paste("Total number of elements:", total_elements))

# Print the proportion of matches
proportion_matches <- num_matches / total_elements
print(paste("Proportion of matches:", proportion_matches))

# compare negative log p
# Rename 'Genes' column in df_sandeep_sorted to 'exon' to match with all_data_sorted
names(df_sandeep_sorted)[names(df_sandeep_sorted) == "Genes"] <- "exon"

# Merge the two dataframes
merged_df <- merge(all_data_sorted, df_sandeep_sorted, by = "exon")

# Calculate -log10(adjusted_pvalues) for both dataframes
merged_df$neg_log_pval_all_data <- -log2(merged_df$adjusted_pvalues.x)
merged_df$neg_log_pval_sandeep <- -log2(merged_df$adjusted_pvalues.y)

# Calculate the abs difference between -log10(adjusted_pvalues) of 'exon' and 'Genes'
merged_df$abs_diff_neg_log_pval <- abs(merged_df$neg_log_pval_all_data - merged_df$neg_log_pval_sandeep)
# Calculate maximum difference
max_diff <- max(merged_df$abs_diff_neg_log_pval, na.rm = TRUE)
merged_df[merged_df$abs_diff_neg_log_pval == max_diff,] # as expected, diff(-logP) is max for the top ranking gene. 

# Plot the distribution of the difference
hist.plot <- ggplot(merged_df, aes(x = abs_diff_neg_log_pval)) +
          geom_histogram(binwidth = 0.1, fill = "blue", color = "white") +  # Adjust binwidth as needed
          geom_vline(aes(xintercept = max_diff), color = "red", linetype = "dashed", linewidth = 1) +
          annotate("text", x = max_diff, y = Inf, label = paste("Max Diff =", round(max_diff, 2)), vjust = 2, hjust = 1.5, color = "red") +
          theme_minimal() +
          labs(x = "Absolute Difference in -log2(bacon P-value)", y = "Count",
               title = paste0("GENESIS vs MMAP (", trait, "): abs diff -log2(bacon P)"))
hist.plot
ggsave(hist.logbaconp.path, plot = hist.plot, width = 10, height = 7, dpi = 300)

# plot (zoom-in version toward the extreme values)
cutoff = 1.5
filtered_df <- merged_df[merged_df$abs_diff_neg_log_pval >= cutoff,]
hist.zoom.plot <- ggplot(filtered_df, aes(x = abs_diff_neg_log_pval)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "white") +  # Adjust binwidth as needed
  geom_vline(aes(xintercept = max_diff), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max_diff, y = Inf, label = paste("Max Diff =", round(max_diff, 2)), vjust = 2, hjust = 1.5, color = "red") +
  theme_minimal() +
  labs(x = paste("Difference in -log2(Adjusted P-value), cutoff=", cutoff), y = "Count",
       title = paste0("GENESIS vs MMAP (", trait, "): zoomed abs diff -log2(bacon P)")) +
  coord_cartesian(xlim = c(1, max(merged_df$abs_diff_neg_log_pval)), ylim = c(0, 20))  # Set x limits and y limits
hist.zoom.plot
# Save the plot
ggsave(hist.logbaconp.zoomed.path, plot = hist.zoom.plot, width = 10, height = 7, dpi = 300)
