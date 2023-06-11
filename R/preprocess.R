#!/usr/bin/env Rscript
library(optparse)
library(dplyr)
library(stringr)
library(diffUTR)
library(edgeR)
library(DESeq2)
library(SummarizedExperiment)
library(foreach)

# Get command-line arguments
parser <- OptionParser()
parser <- add_option(parser, c("--trait"), type="character",
                     default="BMI",
                     help="trait ex) BMI",
                     metavar="trait name")
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

# Set command line input
trait <- opt$trait

# PARAMETERS
visitcode = 1 # 1 or 2

#output file path
out.path <- paste0("outputs/preprocessed/exonCount_", trait, ".RDS")
vst.out.path <- paste0("outputs/preprocessed/vst_", trait, ".RDS")
pc.out.path <- paste0("outputs/preprocessed/pc_", trait, ".csv")

# input file paths
exon.count.path <- "data/diffUTR_count_combined.RDS"
trait.file.path <- paste0("data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_",trait,".csv")
whatdatall.path <- "data/whatdatall.csv"
qualimap.path <- "data/llfs_qualimap_20230323.csv"
gtfPATH = "data/gencode.v38.annotation.gtf"



# read files
phenotype_df <- read.csv(trait.file.path)
whatdatall_df <- read.csv(whatdatall.path)
qualimap_df <- read.csv(qualimap.path)
rds <- as.data.frame(readRDS(exon.count.path))

# Mapping Subject IDs to Sample IDs
mapped_ids <- whatdatall_df %>% 
  semi_join(phenotype_df, by = c("subject" = "subject")) %>% 
  distinct(subject, .keep_all = TRUE) %>% 
  mutate(mapped = if_else(duplicated(id), "Unmapped", "Mapped")) %>% 
  pull(id)

# Filter qualimap df 
qualimap_filtered <- qualimap_df %>% 
  filter(id %in% mapped_ids, str_detect(visit, as.character(visitcode))) %>% 
  arrange(percent_intergenic) %>% 
  distinct(id, .keep_all = TRUE)

# filter for related samples.
filesToKeep = c()
for (bamfile in colnames(rds)){
  # Hard coded based on the filename format of diffUTR runs on 2561 bamfiles
  if (!grepl('pool',bamfile, fixed=T)){ # remove differently formatted bamfiles from the downstream analysis
    sampleID_visitcode <- strsplit(bamfile, split='[.]')[[1]][5] # HARDCODED
    
    sampleID <- str_extract(sampleID_visitcode, '^\\d+')
    
    visitCode <- str_extract(sampleID_visitcode, "visit_\\d|visit\\d|vist\\d|vist_\\d")
    if(grepl(as.character(visitcode), visitCode) & sampleID %in% mapped_ids){ 
      visitCodeToUse = qualimap_filtered[qualimap_filtered$id == sampleID,c('visit')]
      if (visitCode == visitCodeToUse){
        filesToKeep = c(filesToKeep, bamfile)
      }
    }
    if (is.na(visitCode)){
      print(bamfile)
    }
  }
}

# drop unrelavent bamfiles.
rds <- rds[, colnames(rds) %in% filesToKeep]

# rename the columns to s<subjectID>_<visitcode>
subjects = c()
for (bamfile in colnames(rds)){
  sampleID_visitcode <- strsplit(bamfile, split='[.]')[[1]][5]
  
  sampleID <- str_extract(sampleID_visitcode, '^\\d+')
  subjectID <- unique(whatdatall_df[which(whatdatall_df$id == sampleID),]$subject)
  subject <- paste0('s', subjectID, '_v1')
  subjects = c(subjects, subject)
}
# rename the columns
colnames(rds) <- subjects
print(paste("generated RDS file has dimension of ",dim(rds)))
saveRDS(rds, file=out.path)

## filter for bins: QC
bins <- prepareBins(gtfPATH, stranded=FALSE)
exp_se <- SummarizedExperiment(assays=list(counts=rds),
                               rowData=bins)
# Look into bin width and filter
MinBinWidth <- 45 # it was chosen so that the number of the remaining bins < 1M
row_width_filter <- width(ranges(bins)) > MinBinWidth
exp_se <- exp_se[row_width_filter,]

# remove outlying bins. (bins with too many counts)
maxAvgCount = 1000000
myRowMeans <- assays(exp_se)$counts %>% rowMeans(.)
myRowMeans_sorted <- sort(myRowMeans, decreasing=T)
myRowMeans_sorted[1:20]
row_excessive_filter <- myRowMeans < maxAvgCount
exp_se <- exp_se[row_excessive_filter,]

# Filter for lowly expressed exons.
min_num_samples = floor(ncol(exp_se) * 0.015)
expression_level_filter <- rowSums(cpm(assays(exp_se)$counts > 3)) >= min_num_samples
exp_passing <- exp_se[expression_level_filter, ] # no need sample filter as we are using samples used for TWAS

# counts assay is currently a data.frame, need to convert into data.matrix
assays(exp_passing)$counts <- data.matrix(assays(exp_passing)$counts)

# VST transformation
# convert RangedSummarizedExperiment object to DESeqDataSet
dds <- DESeqDataSet(exp_passing, design=~1)
dds <- estimateSizeFactors(dds)
vst_obj <- vst(dds, blind=T)

print(paste0("after vst dim: ", dim(vst_obj)))
# output the
saveRDS(vst_obj, file=vst.out.path)

# get the count data
vst_expr <- assay(vst_obj)
# transpose
vst_expr <- t(vst_expr)

# pc
results <- prcomp(vst_expr, center = TRUE, scale = FALSE)
print("Pc calculated")
print(dim(results$x))

#Write PCs into a matrix
write.csv(results$x, pc.out.path, quote = FALSE)








