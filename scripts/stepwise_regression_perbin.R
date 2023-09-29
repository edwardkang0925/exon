.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript
library(optparse)
library(foreach)
library(stringr)

parser <- OptionParser()
parser <- add_option(parser, c("--combined_df"), type="character",
                     default="outputs/beforeRegression/fhshdl/combined_df_beforeRegression_fhshdl.RDS",
                     help="path to combined_df (covariates matrix) file",
                     metavar="'/path/to/file/combined_df'")
parser <- add_option(parser, c("--exon_expression_path"), type="character",
                     default="outputs/beforeRegression/fhshdl/sliced/exonSlice_9.RDS",
                     help="path to the sliced exon expression RDS file",
                     metavar="'/path/to/dir/containing/exonSlice_<int>.RDS'")
parser <- add_option(parser, c("--output_directory"), type="character",
                     default="outputs/afterRegression/fhshdl/sliced",
                     help="output location",
                     metavar="/path/to/output/dir")


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
residualFilePrefix = "residual_"
# extract slice index by extracting the number between the last underscore and the last dot in the file name
sliceIndex = opt$exon_expression_path %>% str_extract("(?<=_)[^_\\.]+(?=\\.)") %>% as.character()


combined_df <- readRDS(opt$combined_df)
exon_expression_df <- readRDS(opt$exon_expression_path)

residuals_df = data.frame(matrix(ncol = dim(exon_expression_df)[2] , nrow = dim(exon_expression_df)[1]))
rownames(residuals_df) <- rownames(combined_df)
colnames(residuals_df) <- colnames(exon_expression_df)

#check if residuals are normal. Create list of p values and w values for ben shapiro test.
#residuals_normality_shapiro_p_values <- vector(mode='list', length=dim(exon_expression_df)[2])
#residuals_normality_shapiro_w_values <- vector(mode='list', length=dim(exon_expression_df)[2])
base_covariates_list = c("age", "sex", "fc_DK", "fc_NY", "fc_PT", "amon", "aneu", "plt", "rbc", "wbc", "percent_intergenic", "age_2")
additional_covariates_list = colnames(combined_df)[12:42] # plate + pcs


count = 0
index <- 1
temp_residuals <- NULL
#residual adjustment using stepwise regression
# for (values in colnames(exon_expression_df)){
#   #new code start
#   count = count + 1
#   if(count %% 10 == 0){
#     print(count)
#   }
#   combined_df["gene"] <- as.numeric(exon_expression_df[,values])
#   fixed_effects <- as.formula(paste("gene", paste(c(base_covariates_list, additional_covariates_list), collapse="+"), sep="~"))
#   model <- lm(fixed_effects, data=combined_df)
#   tryCatch({
#     model_step <- step(model, scope=list(upper = as.formula(paste("~", paste(c(base_covariates_list, additional_covariates_list), collapse="+"))),
#                                         lower = as.formula(paste("~", paste(base_covariates_list, collapse="+")))), trace=FALSE)
#     temp_residuals <- resid(model_step)
#     residuals_df[,index] <- temp_residuals
#   }, error = function(e) {
#     print(paste("Error in stepwise regression for gene", values))
#     temp_residuals <- runif(dim(residuals_df)[1])
#     residuals_df[,index] <- temp_residuals
#   })
#   #new code end
#   combined_df["gene"] <- NULL
#   #residuals_normality_shapiro_p_values[[index]] <- shapiro.test(temp_residuals)$p.value
#   #residuals_normality_shapiro_w_values[[index]] <- shapiro.test(temp_residuals)$W
#   index <- index + 1
# }
## For test pipeline run on local with sample data, use regression instead of step-wise regression.
for (values in colnames(exon_expression_df)){
  #new code start
  count = count + 1
  if(count %% 10 == 0){
    print(count)
  }
  combined_df["gene"] <- as.numeric(exon_expression_df[,values])
  fixed_effects <- as.formula(paste("gene", paste(c(base_covariates_list), collapse="+"), sep="~"))
  model <- lm(fixed_effects, data=combined_df)

  # Instead of stepwise regression, just compute residuals from the original model
  temp_residuals <- resid(model)
  residuals_df[,index] <- temp_residuals

  #new code end
  combined_df["gene"] <- NULL
  index <- index + 1
}
if (!dir.exists(opt$output_directory)) {
  dir.create(opt$output_directory, recursive = TRUE)
}
print(paste0(opt$output_directory, residualFilePrefix, sliceIndex,".RDS"))
saveRDS(residuals_df, file=paste0(opt$output_directory, residualFilePrefix, sliceIndex,".RDS"))

