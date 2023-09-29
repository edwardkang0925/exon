.libPaths("/project/renv/library/R-4.2/aarch64-unknown-linux-gnu")
#!/usr/bin/env Rscript
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--input_files"), type="character",
                     default="outputs/residuals_slice/",
                     help="list of RDS files to combine",
                     metavar="/path/to/dir/with/RDSfiles")

parser <- add_option(parser, c("--output_prefix"), type="character",
                     default="combined/",
                     help="output location",
                     metavar="/path/to/output/dir")


# main -------------------------------------------------------------------------

opt = parse_args(parser)

filenames <- unlist(strsplit(opt$input_files, ","))
print(filenames)
print(getwd())
# Print to check if the filenames are retrieved correctly
# print(paste("Number of files found:", length(filenames)))

# if(length(filenames) == 0){
#   print(paste("trying to read files from ",opt$input_files))
#   stop("No RDS files found in the specified directory")
# }else{
#   print(filenames)
# }

if (!file.exists(filenames[1])) {
  stop(paste("File does not exist:", filenames[1]))
} else if (file.access(filenames[1], mode = 4) != 0) {
  stop(paste("Do not have read permissions for file:", filenames[1]))
} else {
  combined_obj <- tryCatch({
    readRDS(filenames[1])
  }, error = function(e) {
    stop(paste("Error reading file:", filenames[1], "\n", e))
  })
}

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


