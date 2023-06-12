.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript

# Load necessary library
library(data.table)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")
parser <- add_option(parser, c("--inputDir"), type="character",
                     default="outputs/conditionalGenesisOut_regression_all/",
                     help="path to directory containing RDS files to combine",
                     metavar="'/path/to/dir/with/RDSfiles'")
parser <- add_option(parser, c("--outputFile"), type="character",
                     default="outputs/combinedConditionalGenesis/BMI/comnbinedConditionalGenesis_BMI.RDS",
                     help="output file path",
                     metavar="'/path/to/output/file'")

opt = parse_args(parser)



trait <- opt$trait

# Specify the path to the directory containing the RDS files
path_to_rds_files <- opt$inputDir
outputPATH <- opt$outputFile

# create output dir
directory_path <- dirname(outputPATH)
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}
# Get a list of all RDS files in the directory
rds_files <- list.files(path = path_to_rds_files, pattern = "*.RDS", full.names = TRUE)

# Initialize an empty data table
exons.genesis <- data.table()

# Loop over each RDS file
for (file in rds_files) {
  # Load the data from the RDS file
  data <- readRDS(file)
  print(dim(data))
  # rbind the data to the existing data table
  exons.genesis <- rbindlist(list(exons.genesis, data))
}
# sort the combined result by p-value
exons.genesis <- exons.genesis[order(exons.genesis$p_seq_regression), ]

# Write the combined data table to an RDS file
saveRDS(exons.genesis, file = outputPATH)
