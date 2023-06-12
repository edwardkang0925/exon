.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--input_prefix"), type="character",
                     default="outputs/residuals_slice/",
                     help="location of RDS files to combine",
                     metavar="/path/to/dir/with/RDSfiles")

parser <- add_option(parser, c("--output_prefix"), type="character",
                     default="combined/",
                     help="output location",
                     metavar="/path/to/output/dir")


# main -------------------------------------------------------------------------

opt = parse_args(parser)

filenames <- list.files(opt$input_prefix, pattern="*.RDS", full.names=TRUE)

# Print to check if the filenames are retrieved correctly
# print(paste("Number of files found:", length(filenames)))

# if(length(filenames) == 0){
#   print(paste("trying to read files from ",opt$input_prefix))
#   stop("No RDS files found in the specified directory")
# }else{
#   print(filenames)
# }

combined_obj <- readRDS(filenames[1])

if(length(filenames) > 1){
  for (filename in filenames[c(2:length(filenames))]){
    new_obj <- readRDS(filename)
    combined_obj <- cbind(combined_obj, new_obj)
  }
}

if (!dir.exists(opt$output_prefix)) {
  dir.create(opt$output_prefix, recursive = TRUE)
}

saveRDS(combined_obj, file=paste0(opt$output_prefix, "residualExonExpression.RDS"))


