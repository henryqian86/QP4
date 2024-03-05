library(edgeR)
library(ggplot2)
library(IHW)
library(tidyr)
library(reshape2)
library(DESeq2)

setwd("../data")
df <- read.table("pone.0017820.s009.txt", sep="\t", header = T)

num_genes <- 11529
dispersion <- 0.2 
size <- 1 / dispersion
mean_read_count = 1000
covariate <- rnbinom(n = num_genes, size = size, mu = mean_read_count)
df$covariate = covariate
df = df[c(4, 21)]

calculate_fdr_bh <- function(pvalues, alpha) {
  p_adjusted <- p.adjust(pvalues, method = "BH")
  significant <- p_adjusted < alpha
  estimated_fdr <- mean(p_adjusted[significant])  # Estimate FDR as the mean adjusted p-value for significant tests
  return(estimated_fdr)
}

# subset sizes
subset_sizes <- seq(from = 3000, to = 11000, length.out = 9)

# Nominal alpha levels to evaluate
alpha_levels <- seq(0.01, 0.1, by = 0.01)

# Perform the analysis for different subset sizes
fdr_vs_alpha <- lapply(subset_sizes, function(size) {
  # Sample size rows from the dataframe
  subset_df <- df[sample(nrow(df), size), ]
  
  # Calculate FDR for each alpha level using BH
  fdrs <- sapply(alpha_levels, calculate_fdr_bh, pvalues = subset_df$p.value)
  
  # Return a data frame with the results
  data.frame(SubsetSize = size, Alpha = alpha_levels, FDR_BH = fdrs)
})

# Combine all data frames into one
fdr_vs_alpha_df <- do.call(rbind, fdr_vs_alpha)

# Plot FDR vs. Nominal Alpha for different subset sizes using BH
ggplot(fdr_vs_alpha_df, aes(x = Alpha, y = FDR_BH, color = factor(SubsetSize))) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Add y = x line
  labs(title = "FDR vs. Nominal Alpha(BH)",
       x = "Nominal Alpha Level",
       y = "Estimated FDR",
       color = "Subset Size") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(fdr_vs_alpha_df$Alpha, na.rm = TRUE))) +  # Adjust y-axis if necessary
  scale_x_continuous(limits = c(0, max(fdr_vs_alpha_df$Alpha, na.rm = TRUE)))   # Adjust x-axis if necessary
