# Load necessary library
library(data.table)
library(ggplot2)

# TODO: write slurm script taking trait as an argument or iterate through trait.txt
trait = "BMI"


# Specify the path to the directory containing the RDS files
path_to_rds_files <- paste0("outputs/geneBinCorrelation/slices/")
outputPATH <- paste0("outputs/geneBinCorrelation/combinedCorrelation_",trait,".RDS")
# create output dir
directory_path <- dirname(outputPATH)
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}
# Get a list of all RDS files in the directory
rds_files <- list.files(path = path_to_rds_files, pattern = "*.RDS", full.names = TRUE)

# Initialize an empty data table
corr.df <- data.table()

# Loop over each RDS file
for (file in rds_files) {
  # Load the data from the RDS file
  data <- readRDS(file)
  print(dim(data))
  # rbind the data to the existing data table
  corr.df <- rbindlist(list(corr.df, data))
}
# sort the combined result by p-value
corr.df <- corr.df[order(corr.df$cor), ]

saveRDS(corr.df, file = outputPATH)
# --------------------------------------------- quick plot on correlation -----------------------
corr.hist.plot <- ggplot(corr.df, aes(x = cor)) +
  geom_histogram(aes(y = log10(..count..)), binwidth = 0.01, fill = "blue", color = "white") +  # Adjust binwidth as needed
  theme_minimal() +
  labs(x = "Pearson Correlation (gene-bin)", y = "Log10 Count",
       title = paste0("Correlation between gene and bin: (", trait, ")")) +
  theme(plot.title = element_text(size = 28), 
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 20))

ggsave(paste0("outputs/plots/corr_dist_",trait,".png"), plot = corr.hist.plot, width = 10, height = 7, dpi = 300)


colnames(corr.df)
corr.df[corr.df$cor > 0.95,] # 403 bins

