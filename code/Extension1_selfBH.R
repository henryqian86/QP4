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
  covariate <- effect_sizes + rnorm(num_tests, sd = 0.5)
  
  # Return a dataframe
  data <- data.frame(p.value = pvalues, covariate = covariate, is_alt = is_alt)
  return(data)
}

# Simulate data
num_tests <- 10000
df <- as.data.frame(simulate_informative_data(num_tests))

# Define the sizes of subsets you want to analyze
subset_sizes <- seq(from = 3000, to = 10000, length.out = 8)

# Nominal alpha levels to evaluate
alpha_levels <- seq(0.01, 0.1, by = 0.01)

# Perform the analysis for different subset sizes
fdr_vs_alpha <- lapply(subset_sizes, function(size) {
  # Sample 'size' rows from the dataframe
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
