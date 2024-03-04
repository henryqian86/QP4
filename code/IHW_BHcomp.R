library(edgeR)
library(ggplot2)
library(IHW)
library(tidyr)
library(reshape2)
library(DESeq2)

setwd("/Users/henryqian/Desktop/QP4/data")
df <- read.table("pone.0017820.s009.txt", sep="\t", header = T)

set.seed(1109)
num_genes <- 11529

# Parameters for the negative binomial distribution
dispersion <- 0.2       # Dispersion parameter, where variance = mean + mean^2 * dispersion

# Size parameter for rnbinom (related to the variance)
size <- 1 / dispersion

alphaDis <- function(size, mean_read_count) {
  # Generate simulated read counts
  cov <- rnbinom(n = num_genes, size = size, mu = mean_read_count)
  
  #cov <- df$Low_read_count
  pvalues <- df$p.value
  
  # Running IHW
  alphas <- seq(0.05, 0.10, by = 0.01)
  discoveries <- sapply(alphas, function(alpha) {
    ihw_result <- ihw(pvalues, cov, alpha)
    list(
      IHW = sum(ihw_result@df$adj_pvalue < alpha),
      BH = sum(p.adjust(pvalues, method = "BH") < alpha)
    )
  })
  discoveries_df <- as.data.frame(t(discoveries))
  discoveries_df <- data.frame(lapply(discoveries_df, as.integer))
  
  
  discoveries_df1 <- data.frame(CombinedColumn = c(discoveries_df$BH, discoveries_df$IHW))
  discoveries_df1$Alpha <- rep(alphas,2)
  discoveries_df1$Method <- rep(c("BH", "IHW"), each=length(alphas))
  
  
  # Plotting
  ggplot(discoveries_df1, aes(x = Alpha, y = CombinedColumn, color = Method)) +
    geom_line() + # Add lines
    geom_point() + # Add points
    theme_minimal() + # Use a minimal theme
    labs(x = "Alpha", y = "Discoveries", title = "Alpha vs. Discoveries", color = "Type") +
    scale_color_manual(values = c("IHW" = "blue", "BH" = "red")) # Custom colors
}

alphaDis(1/0.2, 1000)
