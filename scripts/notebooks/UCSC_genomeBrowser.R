library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(tidyr)  # For the gather function used later
library(biomaRt)
# read GRange object containing bins of our interest
df = readRDS("data/pipelineResult/comnbinedConditionalGenesis_lnTG.RDS")

# read bin coordinate
bins <- readRDS("data/bins.RDS")

# subset for conditionally significant bins.
df_sig <- df[df$p_lrt <= 0.01 / dim(df)[1],]

# now subset bins based on bin's name
df_list_sigbins <- list()
for(bin in df_sig$exon){
  # Convert bin name to gene name by removing the last ".X" part
  gene_name <- unlist(strsplit(bin, "\\."))[1:(length(unlist(strsplit(bin, "\\."))) - 1)] %>% paste(collapse = ".")
  bin_number <- as.numeric(tail(unlist(strsplit(bin, "\\.")), 1))
  print(gene_name)
  print(bin_number)
  # Filter bins and store the filtered dataframe in df_list_sigbins
  df_list_sigbins[[length(df_list_sigbins) + 1]] <- bins[bins$gene == gene_name & bins$bin_id == bin_number,]
}
targets <- do.call(c, df_list_sigbins)

# REMOVE later: currently, bins overlapping with multiple genes are problematic as their condtional analysis was not conducted correctly.
# drop bins that have ambiguous gene
targets <- targets[targets$geneAmbiguous == FALSE,]
targets$gene_name <- as.character(targets$gene_name)

cma_tg_genes <- c("AKAP12", "CA8", "CPA3", "ENPP3", "FCER1A", "GATA2", "GCSAML", "HDC", "HRH4", "LINC02458", "MS4A2", "MS4A3", "NTRK1", "SLC45A3")
intersecting_genes <- c("MS4A2", "NTRK1", "AKAP12", "CA8", "GATA2", "HDC", "SLC45A3", "ENPP3")
targets <- targets[targets$gene_name %in% intersecting_genes]
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


