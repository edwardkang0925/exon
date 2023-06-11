library(ggplot2)
# calculate GIF
inflation <- function(p) {
  chisq <- qchisq(1 - p, df = 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

# outname: output path + name for the plot png
# main name: plot title
# pvals: p values as a vector
plot_qq <- function(outname, main_name, pvals, threshold) {
  cex_point <- 1
  observed <- sort(pvals)
  lobs <- -(log10(observed))
  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected) + 1)))
  
  png(outname)
  plot(lexp, lobs, pch = 0,
       cex = cex_point,
       xlab = expression(Expected ~ ~-log[10](italic(p))),
       ylab = expression(Observed ~ ~-log[10](italic(p))),
       font.lab = 2,
       cex.lab = 1,
       cex.axis = 1,
       font.axis = 2,
       main = main_name)
  
  abline(0, 1, col = "red", lwd = 1)
  
  # Add horizontal line at specified threshold
  abline(h = -(log10(threshold)), col = "blue", lwd = 1)
  
  gif <- inflation(pvals)
  legend("topleft",
         legend = c(paste0("GIF = ", format(round(gif, 3), nsmall = 2))),
         ncol = 1,
         bty = "o",
         box.lwd = 1,
         pch = 0,
         cex = 1,
         text.font = 2)
  
  dev.off()
}

plot_histogram <- function(trait, df, column_name, bin_width, output_path) {
  # Check if the column exists in the dataframe
  if (!column_name %in% colnames(df)) {
    stop(paste("Column", column_name, "does not exist in the dataframe."))
  }
  
  # Check if the column contains numeric data
  if (!is.numeric(df[[column_name]])) {
    stop(paste("Column", column_name, "is not numeric."))
  }
  
  # Check if the bin_width is a positive number
  if (!is.numeric(bin_width) || bin_width <= 0) {
    stop("bin_width must be a positive number.")
  }
  
  # Create the histogram
  p <- ggplot(df, aes_string(column_name)) +
    geom_histogram(aes(y = log10(..count..)), binwidth = bin_width, fill = "blue", color = "black") +
    xlab(column_name) +
    ylab("Log Frequency") +
    ggtitle(paste(column_name, ": Checking Uniformity ", trait)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 28), 
          axis.title = element_text(size = 24), 
          axis.text = element_text(size = 20))
  
  # Save the plot to the specified path
  ggsave(filename = output_path, plot = p)
}



traits <- c("BMI", "fhshdl")
for(trait in traits){
  conditional.genesis.result.path <- paste0("outputs/combinedConditionalGenesis/",trait,"/combinedConditionalGenesis_",trait,".RDS")
  conditional.exons.genesis <- readRDS(conditional.genesis.result.path)
  plot_qq(paste0("outputs/plots/conditionalAll/", trait,"_qq.png"),
          paste0(trait," qq plot LRT"),
          conditional.exons.genesis$p_lrt,
          threshold=(0.01 / dim(conditional.exons.genesis)[1]))
  plot_histogram(trait, conditional.exons.genesis, "p_lrt",0.01, paste0("outputs/plots/conditionalAll/", trait, "_pval_hist.png"))
}

for(trait in traits){
  conditional.genesis.result.path <- paste0("outputs/combinedConditionalGenesis/",trait,"/combinedConditionalGenesis_",trait,"_regression.RDS")
  conditional.exons.genesis <- readRDS(conditional.genesis.result.path)
  plot_qq(paste0("outputs/plots/conditionalAll/", trait,"_regression_qq.png"),
          paste0(trait," qq plot SR"),
          conditional.exons.genesis$p_seq_regression,
          threshold=(0.01 / dim(conditional.exons.genesis)[1]))
  plot_histogram(trait, conditional.exons.genesis, "p_seq_regression" ,0.01, paste0("outputs/plots/conditionalAll/", trait, "_sr_pval_hist.png"))
}



