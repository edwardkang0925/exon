library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
# read GRange object containing bins of our interest
targets = readRDS("data/postAnalysis_bins_Grange_objs/bmi_condSigBins.RDS")

# Set up biomaRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Drop unused factor levels from the 'gene' column
mcols(targets)$gene <- factor(mcols(targets)$gene, levels = unique(mcols(targets)$gene))

# Split targets by gene
gene_splits <- split(targets, mcols(targets)$gene)
bins_index = 1
# Loop through each gene subset
for (gene in names(gene_splits)) {
  # Start browser session
  session <- browserSession("UCSC")

  gene_bins <- gene_splits[[gene]]
  gene_name_no_version <- sub("\\..*", "", gene)

  # Fetch gene's coordinates
  gene_info <- getBM(attributes = c("chromosome_name", "start_position", "end_position"),
                     filters = 'ensembl_gene_id',
                     values = gene_name_no_version,
                     mart = mart)

  # Define the range for the entire gene
  gene_range <- GRanges(seqnames = paste0("chr", gene_info$chromosome_name),
                        ranges = IRanges(start = gene_info$start_position,
                                         end = gene_info$end_position))

  targetRanges <- ranges(gene_bins)  # <- Important: Use 'gene_bins' here, not 'targets'

  targetTrack <- with(gene_bins,     # <- Use 'gene_bins' here as well
                      GRangesForUCSCGenome("hg38", chrom = seqnames, ranges = targetRanges))

  # Add bins to UCSC session
  track_name <- paste0("sigBins_", bins_index)
  bins_index = bins_index + 1
  track(session, track_name) <- targetTrack
  browserView(session, targetTrack * -50)

}


