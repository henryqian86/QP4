library(edgeR)
library(ggplot2)
library(IHW)

setwd("../data")
df <- read.table("pone.0017820.s009.txt", sep="\t", header = T)

set.seed(123)
num_genes <- 11529
dispersion <- 0.2 
size <- 1 / dispersion
mean_read_count = 1000
covariate <- rnbinom(n = num_genes, size = size, mu = mean_read_count)
df$covariate = covariate
df = df[c(4, 21)]

# Function to apply both IHW and BH on a subset of data and calculate discoveries
apply_methods <- function(subset_df) {
  pvalues = subset_df$p.value
  cov = subset_df$covariate
  
  # Apply IHW
  ihw_result <- ihw(pvalues, cov, alpha = 0.01)
  discoveries_ihw <- sum(ihw_result@df$adj_pvalue <= 0.01)
  
  # Apply BH
  pvalues_adj_bh <- p.adjust(pvalues, method = "BH")
  discoveries_bh <- sum(pvalues_adj_bh <= 0.01)
  
  # Return the discoveries for both methods
  return(list(discoveries_ihw = discoveries_ihw, discoveries_bh = discoveries_bh))
}

# Function to create subsets and apply both methods
analyze_subsets <- function(df, subset_sizes) {
  results <- list()
  
  for (size in subset_sizes) {
    # Sample size rows from the dataframe
    subset_df <- df[sample(nrow(df), size), ]
    
    # Apply IHW and BH to the subset
    methods_result <- apply_methods(subset_df)
    
    # Store the result
    results[[paste("Size", size)]] <- c(methods_result, SubsetSize = size)
  }
  
  return(results)
}

subset_sizes <- seq(from = 3000, to = 11000, length.out = 9)

# Perform the analysis
results <- analyze_subsets(df, subset_sizes)

# Convert results to a data frame
results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)

# Convert character columns to numeric
results_df$SubsetSize <- as.numeric(as.character(results_df$SubsetSize))
results_df$discoveries_ihw <- as.numeric(as.character(results_df$discoveries_ihw))
results_df$discoveries_bh <- as.numeric(as.character(results_df$discoveries_bh))

# Plot Number of Discoveries vs. Subset Size for both IHW and BH
ggplot() +
  geom_line(data = results_df, aes(x = SubsetSize, y = discoveries_ihw, colour = "IHW")) +
  geom_point(data = results_df, aes(x = SubsetSize, y = discoveries_ihw, colour = "IHW")) +
  geom_line(data = results_df, aes(x = SubsetSize, y = discoveries_bh, colour = "BH")) +
  geom_point(data = results_df, aes(x = SubsetSize, y = discoveries_bh, colour = "BH")) +
  labs(title = "Number of Discoveries vs. Subset Size with RNA-seq dataset",
       x = "Subset Size",
       y = "Number of Discoveries",
       color = "Method") +
  theme_minimal()
