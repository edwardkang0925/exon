library(tidyverse)

# FIXME: trait
trait <- "BMI"
traits <- c("BMI", 'fhshdl') # traits with results
logPath <- "outputs/postConditionalGenesis_all.txt"

for(trait in traits){
  
  # input paths
  conditional.genesis.result.path <- paste0("outputs/combinedConditionalGenesis/",trait,"/combinedConditionalGenesis_",trait,"_regression.RDS")
  
  # load data
  conditional.genesis.df <- readRDS(conditional.genesis.result.path)
  conditional.genesis.df$nlog10p <- -log10(conditional.genesis.df$p_seq_regression)
  nlog10p.max <- max(conditional.genesis.df$nlog10p)
  
  # using two thresholds to separate super sig bins and 
  df.nominal.sig <- conditional.genesis.df[conditional.genesis.df$p_seq_regression < 0.01 / dim(conditional.genesis.df)[1],]
  
  # output # sigbins using two different thresholds.
  write(paste("[",trait,"]",
              dim(df.nominal.sig)[1], "sig bins (bin-wide Bonferroni FDR = 0.01)"), logPath, append=TRUE)
  
  
  ## Make -log10P scatter plot
  # Plot the distribution of the -log10P (per trait)
  nlogp.hist.plot <- ggplot(conditional.genesis.df, aes(x = nlog10p)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "white") +  # Adjust binwidth as needed
    geom_vline(aes(xintercept = nlog10p.max), color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = nlog10p.max, y = Inf, label = paste("Max =", round(nlog10p.max, 2)), vjust = 2, hjust = 1.5, color = "red") +
    theme_minimal() +
    labs(x = "-log10(p_seq_regression)", y = "Count",
         title = paste0("Conditional Genesis all: (", trait, " SR) -log10P")) +
    theme(plot.title = element_text(size = 28), 
          axis.title = element_text(size = 24), 
          axis.text = element_text(size = 20))
  
  ggsave(paste0("outputs/plots/conditionalAll/hist_nlogp_conditionalGenesis_",trait,"_sr_all.png"), plot = nlogp.hist.plot, width = 10, height = 7, dpi = 300)
  
  # nominal sig only: Plot the distribution of the -log10P (per trait)
  nlogp.sig.hist.plot <- ggplot(df.nominal.sig, aes(x = nlog10p)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "white") +  # Adjust binwidth as needed
    geom_vline(aes(xintercept = nlog10p.max), color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = nlog10p.max, y = Inf, label = paste("Max =", round(nlog10p.max, 2)), vjust = 2, hjust = 1.5, color = "red") +
    theme_minimal() +
    labs(x = "-log10(p_seq_regression), zoomed", y = "Count",
         title = paste0("Conditional Genesis all: (", trait, " SR) -log10P")) +
    theme(plot.title = element_text(size = 28), 
          axis.title = element_text(size = 24), 
          axis.text = element_text(size = 20))
  
  ggsave(paste0("outputs/plots/conditionalAll/hist_nlogp_conditionalGenesis_",trait,"_sr_sig_all.png"), plot = nlogp.sig.hist.plot, width = 10, height = 7, dpi = 300)
  
}
