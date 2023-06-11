
rds1 <- readRDS("outputs/afterRegression/ABI/combined/residualExonExpression.RDS")
dim(rds1)
not_ENSG_cols <- grep("ENSG", colnames(rds1), invert = TRUE, value = TRUE)
not_ENSG_cols

rds2 <- readRDS("outputs/genesisInput/ABI/exon_data_001.RDS")
"ABI" %in% colnames(rds2)
rownames(rds2)


kinship.matrix.full.path <- "data/kins.Rdata"
load(kinship.matrix.full.path) # loaded as 'kins'
saveRDS(kins, "data/kinship.full.RDS")

rds3 <- readRDS("outputs/genesisOut/fhshdl/results_027.RDS")
head(rds3)
dim(rds3)


mmap.fvc.df <- read.csv("data/MMAP_INPUT_FINAL_ROUND/twas_input_visit_1_adjusted_fhshdl.csv")
dim(mmap.fvc.df)
"ENSG00000179348.12" %in% colnames(mmap.fvc.df)
colnames(mmap.fvc.df)

# calculate GIF of GENESIS gene level on fhshdl 
# calculate GIF
inflation <- function(p) {
  chisq <- qchisq(1 - p, df = 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
df <- read.csv("expected/combined/genesisCombinedGeneLevel.csv")
dim(df)
colnames(df)
GIF_genesis_fhshdl <- inflation(df$pval)
