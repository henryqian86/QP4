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


summary_stats <- lapply(names(results), function(size) {
  data <- results[[size]]@df
  discoveries <- sum(data$adj_pvalue <= 0.05)
  return(data.frame(SubsetSize = as.numeric(sub("Size ", "", size)), 
                    Discoveries = discoveries))
})

# Combine all summaries into a single dataframe
summary_df <- do.call(rbind, summary_stats)

# Plot Number of Discoveries vs. Subset Size
ggplot(summary_df, aes(x = SubsetSize, y = Discoveries)) +
  geom_line() +
  geom_point() +
  labs(title = "Number of Discoveries vs. Subset Size",
       x = "Subset Size",
       y = "Number of Discoveries") +
  theme_minimal()


