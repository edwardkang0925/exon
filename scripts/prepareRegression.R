.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript
library(optparse)
library(foreach)
library(MASS)
library(here)
library(tidyr)
library(dplyr)
library(stringr)
library(fastDummies)
# library(vsn)
library(hexbin)
library(ggplot2)
library(reshape2)
library(DESeq2)

# Get command-line arguments
parser <- OptionParser()
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")
parser <- add_option(parser, c("--number_of_plates"), type="numeric", default=28,
                      help="number of plates")
parser <- add_option(parser, c("--visitcode"), type="character", default="v1",
                      help="visitcode ex) v1 or v2")
parser <- add_option(parser, c("--numSlice"), type="numeric", default=100,
                      help="number of slices")
parser <- add_option(parser, c("--pcFilePath"), type="character", default="",
                      help="pc file path")
parser <- add_option(parser, c("--vstFilePath"), type="character", default="",
                      help="vst file path")
opt = parse_args(parser)

# check cmd line arguments
x = foreach(
  i = names(opt)
) %do% {
  input_value=opt[[i]]
  if(input_value==''){
    stop(sprintf("ARGUMENT --%s IS REQUIRED",i))
  }else if(i %in% c('dna','rna')){
    if(!file.exists(input_value)){
      stop(sprintf("FILE %s DOES NOT EXIST",input_value))
    }
  }
}

# after preprocess_appending.R (subset samples for trait, bin filter, vst, pc)
# command line argument
phenotype= opt$trait
number_of_plates= opt$number_of_plates
visitcode_param = opt$visitcode # "v1"
numSlice = opt$numSlice # 100
# create output path
outputPATH = paste0("outputs/beforeRegression/",phenotype,"/")
slicedOutputPATH = paste0(outputPATH, "sliced/")

# create output path directory if not existing.
if (!dir.exists(slicedOutputPATH)) {
  dir.create(slicedOutputPATH, recursive = TRUE)
}

# input files paths
misc_file_path="/project/data/MISC_FINAL_ROUND/"
input_file_path="outputs/preprocessed/" #FIXME: when using nextflow, take this as command line argument.
phenotype_file_path="/project/data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/"

blood_all_filepath <- paste0(misc_file_path, "blood_all.csv")
metadata_filepath <- paste0(misc_file_path, "sample_meta_20220608.csv")
subject_age_fc_filepath <- paste0(misc_file_path, "subject_fc_sex_age_revised_no_blood_cancer.csv")
phenotype_filepath <- paste0(phenotype_file_path, "adjusted_",phenotype,".csv")
pcs_filepath <- paste0(opt$pcFilePath)
exon_expression_filepath <- paste0(opt$vstFilePath)

#reading and pre-processing_files
blood_all_df <- read.table(file = blood_all_filepath,  header = TRUE, sep = ",")
metadata_df <- read.table(file = metadata_filepath, header = TRUE, sep = ",")
subject_age_fc_df <- read.table(file = subject_age_fc_filepath, header = TRUE, sep = ",")
phenotype_df <- read.table(file = phenotype_filepath, header = TRUE, sep = ",")
exon_expression_df <- readRDS(file = exon_expression_filepath)
print(dim(exon_expression_df))
pcs_df <- read.table(file = pcs_filepath, header = TRUE, sep = ",")
print(dim(pcs_df))
#remove a row with fc null input
subject_age_fc_df <- subject_age_fc_df[subject_age_fc_df$fc != "", ]

# create covariates vector
covariates.list <- c("subject", "percent_intergenic")
for (i in c(2:number_of_plates)){
  covariates.list <- c(covariates.list, paste0("plate_", as.character(i)))
}

## preprocess pcs_df 
pcs.10.df <- pcs_df[, c("X", "PC1", "PC2", "PC3", "PC4")]
pcs.10.df <- dplyr::rename(pcs.10.df, subject = X) # rename column
# reformat subject column: s<subjectID>_v<visitcode> -> <subjectID>
pcs.10.df$subject <- str_extract(pcs.10.df$subject, "\\d+")


## preprocess metadata_df (subsetting, reformatting)
metadata_df <- dummy_cols(metadata_df, select_columns = "plate", remove_first_dummy = TRUE, remove_selected_columns = TRUE)
metadata_df <- metadata_df %>%
  filter(str_detect(subject_count_headers, visitcode_param)) %>%
  filter(subject %in% pcs.10.df$subject) %>%
  dplyr::select(any_of(covariates.list)) 


## preprocess subject_age_fc_df
subject_age_fc_df <- dummy_cols(subject_age_fc_df, select_columns = "fc", remove_first_dummy = TRUE, remove_selected_columns = TRUE)
subject_age_fc_df <- subject_age_fc_df %>%
  filter(visitcode == as.integer(str_extract(visitcode_param, "\\d+"))) %>%
  filter(subject %in% pcs.10.df$subject)
dim(subject_age_fc_df)


## preprocess blood_all_df
blood_covariates <- c("subject", "amon", "aneu", "rbc", "wbc", "plt")
blood_all_df <- blood_all_df %>%
  filter(str_detect(visitcode, str_extract(visitcode_param, "\\d+"))) %>%
  filter(subject %in% pcs.10.df$subject) %>%
  dplyr::select(any_of(blood_covariates)) 
dim(blood_all_df)

## combine preprocessed dataframes
dfs.to.combine <- list(subject_age_fc_df,blood_all_df, metadata_df, pcs.10.df)
# Function to change "subject" column to numeric
change_subject_to_numeric <- function(df) {
  df[["subject"]] <- as.numeric(df[["subject"]])
  return(df)
}
# change datatype of subject to numeric. Some were character.
dfs.to.combine <- lapply(dfs.to.combine, change_subject_to_numeric)
# merge dataframes by subject
combined_df <- dfs.to.combine %>% purrr::reduce(full_join, by = "subject")
dim(combined_df)
# replace empty string as NA
combined_df[combined_df==""]<-NA
combined_df <- combined_df %>% drop_na()
dim(combined_df)
# change blood count to numeric, add rownames, sort by subjectID, add age squared,
combined_df <- combined_df %>%
  mutate_at(vars(plt, rbc, wbc), ~as.numeric(as.character(.))) %>%
  tibble::column_to_rownames("subject") %>%
  arrange(as.numeric(row.names(.))) %>%
  mutate(age_2 = age ^ 2) %>%
  dplyr::select(-visitcode)
dim(combined_df)

## Modify the exon_expression_df: transpose, modify subject id to numeric.
# Convert the DESeqTransform to a matrix
matrix_data <- assay(exon_expression_df)
# Create a dataframe
exon_expression_df2 <- as.data.frame(matrix_data)
# Transpose, rename subjectID, sort by subjectID, and set rownames == subjectID
transposed_df <- exon_expression_df2 %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rowname") %>% # Move the rownames to a column
  mutate(rowname = as.numeric(str_extract(rowname, "\\d+(?=_v)"))) %>%
  arrange(rowname) %>% # Sort by subject ID
  tibble::column_to_rownames("rowname") # Move the subject ID back to rownames

saveRDS(combined_df, file=paste0(outputPATH, "combined_df_beforeRegression_",phenotype,".RDS"))
saveRDS(phenotype_df, file=paste0(outputPATH, "phenotype_df_beforeRegression_",phenotype,".RDS"))

## Output sliced exon_expression for array jobs. Stepwise regression takes long time
# Identify the column names that are exon ids
exon_cols <- grep("ENSG", colnames(transposed_df), value = TRUE)
# Calculate the number of exons per file
exons_per_file <- ceiling(length(exon_cols) / numSlice)
# Split the exons into numSlice groups
split_exons <- split(exon_cols, ceiling(seq_along(exon_cols) / exons_per_file))
# Create each sliced RDS file
for(i in seq_along(split_exons)) {
  file_data <- transposed_df[, split_exons[[i]], drop = FALSE]
  saveRDS(file_data, file=paste0(slicedOutputPATH,"exonSlice_", i,".RDS"))
}

