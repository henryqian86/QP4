library(ggplot2)
library(IHW)

set.seed(42)  # Setting seed for reproducibility

simulate_informative_data <- function(num_tests, prop_alt = 0.2, effect_size = 2) {
  # Simulate effect sizes
  effect_sizes <- rnorm(num_tests, sd = 1)
  # Determine which tests are alternatives based on proportion
  is_alt <- runif(num_tests) < prop_alt
  # Add an effect to the alternatives
  effect_sizes[is_alt] <- effect_sizes[is_alt] + effect_size
  
  # Simulate p-values from effect sizes
  # Assuming a two-sided test, and under the null, effect sizes follow a standard normal distribution
  pvalues <- 2 * pnorm(-abs(effect_sizes))
  
  # Create an informative covariate
  # Here, the covariate is correlated with effect size, with some noise
  covariate <- runif(n=num_tests,min= -4,max=7)
  
  # Return a dataframe
  data <- data.frame(p.value = pvalues, covariate = covariate, is_alt = is_alt)
  return(data)
}

# Simulate data
num_tests <- 10000
df <- as.data.frame(simulate_informative_data(num_tests))

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

subset_sizes <- seq(from = 3000, to = 10000, length.out = 8)

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
