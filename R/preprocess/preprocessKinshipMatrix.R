traits <- readLines("trait.txt")
for (trait in traits){
  exon.path <- paste0("outputs/genesisInput/",trait,"/exon_data_001.RDS")
  output.path <- "data/kinship/"
  
  # path to required files: 
  kinship.matrix.full.path <- "data/kinship.full.RDS"
  kinship_matrix <- readRDS(kinship.matrix.full.path)
  
  # Load the exon expression data
  exon_data <- readRDS(exon.path)
  not_ENSG_cols <- grep("ENSG", colnames(exon_data), invert = TRUE, value = TRUE)
  not_ENSG_cols
  
  # subset kinship matrix to match the order of subjects with the exon_data: 
  # * ASSUMPTION: the order of subjects across exon_data_001 to exon_data_100 is consistent because we are using 001 as a reference. 
  common_ids <- intersect(rownames(exon_data), rownames(kinship_matrix))
  kinship_matrix <- kinship_matrix[common_ids, common_ids]
  kinship_matrix <- kinship_matrix[rownames(exon_data), rownames(exon_data)]
  
  saveRDS(kinship_matrix, paste0(output.path, "kinship_", trait, ".RDS"))
}