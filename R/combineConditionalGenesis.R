# Load necessary library
library(data.table)

# TODO: write slurm script taking trait as an argument or iterate through trait.txt
trait = "BMI"
traits = c("BMI", "fhshdl")
for(trait in traits){
  # Specify the path to the directory containing the RDS files
  path_to_rds_files <- paste0("outputs/conditionalGenesisOut_regression_all/",trait)
  outputPATH <- paste0("outputs/combinedConditionalGenesis/",trait,"/combinedConditionalGenesis_",trait,"_regression.RDS")
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
  
  # # Bonferroni correction on nominal p
  # exons.genesis$p_value_bonferroni <- p.adjust(exons.genesis$pval, method = "bonferroni")
  # 
  # # BACON correction
  # genesis_z_score = qnorm(exons.genesis$pval/2, lower.tail=FALSE) *
  #   sign(exons.genesis$beta / exons.genesis$se)
  # genesis_bacon <- bacon(teststatistics = genesis_z_score, niter=10000L, nburnin = 3000L)
  # exons.genesis$pval_bacon <- pval(genesis_bacon, corrected = TRUE)
  # print(dim(exons.genesis))
  # save the combined result.
  saveRDS(exons.genesis, file = outputPATH)
}
