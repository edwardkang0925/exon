library(data.table)

# split_exon_files.R
exon.path <- "expected/input/twas_input_visit_1_adjusted_fhshdl.csv"

# Load the exon expression data
exon_data <- data.table::fread(exon.path)

# Identify the column names that are exon ids
exon_cols <- grep("ENSG", colnames(exon_data), value = TRUE)

# Calculate the number of exons per file
exons_per_file <- ceiling(length(exon_cols) / 100)
print(exons_per_file)
# Split the exons into 200 groups
split_exons <- split(exon_cols, ceiling(seq_along(exon_cols) / exons_per_file))

# Create each file
for(i in seq_along(split_exons)) {
  file_data <- exon_data[, c("EGO","fhshdl", split_exons[[i]]), with = FALSE]
  data.table::fwrite(file_data, sprintf("expected/splitted/fhshdl/exon_data_%03d.csv", i))
}
