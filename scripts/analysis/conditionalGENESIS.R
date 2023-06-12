library(GENESIS)

#trait <- "lnTG"
traits <- c("lnTG", "adjTC", "BMI", "fhshdl") # traits with sig genes

for(trait in traits){
  # input path (generated from preprocessForComparison.R)
  input.dir <- paste0("outputs/conditionalGenesisInput/", trait)
  kinship.matrix.path <- paste0("data/kinship/kinship_",trait,".RDS")
  trait.df.path <- paste0("data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_", trait,".csv")
  output.path <- paste0("outputs/conditionalGenesisOut/ConditionalGenesis_",trait,".RDS")
  
  # read input data
  kinship_matrix <- readRDS(kinship.matrix.path)
  trait.df <- read.csv(trait.df.path)
  
  # Get the list of RDS files
  rds_files <- list.files(input.dir, pattern = "\\.RDS$", full.names = TRUE)
  
  # initialize per trait summary file
  results <- data.frame()
  
  # Loop through each RDS file (samples by 1 sig gene, bins within the gene)
  for (rds_file in rds_files) {
    # Load the RDS file
    data <- readRDS(rds_file)
    # append trait FIXME: could append it when preprocessing.
    data$subject <- rownames(data)
    data.subject.trait.appended <- merge(data, trait.df, by="subject") %>%
      column_to_rownames(var="subject")
  
    # Find columns to work with
    colnames_data <- colnames(data.subject.trait.appended)
    gene_column <- grep("^ENSG\\d+\\.\\d+$", colnames_data, value = TRUE)
    exon_columns <- grep("^ENSG\\d+\\.\\d+\\.\\d+$",colnames_data, value = TRUE)
    
    
    # Loop through each exon: do trait ~ exon | gene considering kinship_matrix
    for (exon_col in exon_columns) {
      model_data <- data.frame(trait = data.subject.trait.appended[[trait]],  
                               exon = data.subject.trait.appended[[exon_col]],
                               gene = data.subject.trait.appended[[gene_column]],
                               row.names = rownames(data.subject.trait.appended))
      # Fit the null model
      null_model <- fitNullModel(
        x = model_data,
        outcome = "trait",
        covars = "gene",  
        cov.mat = kinship_matrix
      )
      
      # Fit the alternative model and perform the association test
      alt_model <- fitNullModel(
        x = model_data,
        outcome = "trait",
        covars = c("gene","exon"), 
        cov.mat = kinship_matrix
      )
      
      # calculate LRT
      lrt <- 2*(alt_model$logLik- null_model$logLik)
      p_value_lrt <- pchisq(2 * lrt, df = 1, lower.tail = FALSE)
      p_value_lrt_bon <- min(1.0, p_value_lrt * length(exon_columns))
      
      # Append the results to the results data frame
      results <- rbind(results, data.frame(trait = trait, exon = exon_col, p_lrt = p_value_lrt, p_lrt_bon = p_value_lrt_bon))
    }
  }
  # save results
  saveRDS(results, file = output.path)
}  


