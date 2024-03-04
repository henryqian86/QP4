library(edgeR)
library(ggplot2)
library(IHW)
library(tidyr)
library(reshape2)
library(DESeq2)

setwd("/Users/henryqian/Desktop/QP4/data")
df <- read.table("pone.0017820.s009.txt", sep="\t", header = T)

num_genes <- 11529
dispersion <- 0.2 
size <- 1 / dispersion
mean_read_count = 1000
covariate <- rnbinom(n = num_genes, size = size, mu = mean_read_count)
df$covariate = covariate
df = df[c(4, 21)]


# Function to apply IHW on a subset of data
apply_ihw <- function(subset_df) {
  pvalues = subset_df$p.value
  cov = subset_df$covariate
  ihw_result <- ihw(pvalues, cov, alpha = 0.05)
  return(ihw_result)
}

# Function to create subsets and apply IHW
analyze_subsets <- function(df, subset_sizes) {
  results <- list()
  
  for (size in subset_sizes) {
    # Sample `size` rows from the dataframe
    subset_df <- df[sample(nrow(df), size), ]
    
    # Apply IHW to the subset
    ihw_result <- apply_ihw(subset_df)
    
    # Store the result
    results[[paste("Size", size)]] <- ihw_result
  }
  
  return(results)
}

subset_sizes <- seq(from = 3000, to = 11000, length.out = 9)

# Perform the analysis
results <- analyze_subsets(df, subset_sizes)


# Function to calculate FDR for a given alpha level
calculate_fdr <- function(adj_pvalues, alpha) {
  significant <- adj_pvalues < alpha
  estimated_fdr <- mean(adj_pvalues[significant])  # Estimate FDR as the mean adjusted p-value for significant tests
  return(estimated_fdr)
}

# Nominal alpha levels to evaluate
alpha_levels <- seq(0.01, 0.1, by = 0.01)

# Assuming 'results' is your list of ihw_result objects
fdr_vs_alpha <- lapply(names(results), function(size) {
  adj_pvalues <- results[[size]]@df$adj_pvalue
  fdrs <- sapply(alpha_levels, calculate_fdr, adj_pvalues = adj_pvalues)
  data.frame(Alpha = alpha_levels, FDR = fdrs, SubsetSize = size)
})

# Combine all data frames into one
fdr_vs_alpha_df <- do.call(rbind, fdr_vs_alpha)

# Plot FDR vs. Nominal Alpha for different subset sizes
ggplot(fdr_vs_alpha_df, aes(x = Alpha, y = FDR, color = SubsetSize)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Add y = x line
  labs(title = "FDR vs. Nominal Alpha (IHW)",
       x = "Nominal Alpha Level",
       y = "Estimated FDR",
       color = "Subset Size") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(fdr_vs_alpha_df$Alpha, na.rm = TRUE))) +  # Adjust y-axis if necessary
  scale_x_continuous(limits = c(0, max(fdr_vs_alpha_df$Alpha, na.rm = TRUE)))   # Adjust x-axis if necessary