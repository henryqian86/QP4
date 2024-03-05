library(ggplot2)
library(IHW)
library(FactoMineR)

set.seed(42)

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
  covariate2 <- 2*effect_sizes + rnorm(num_tests, sd = 0.7)
  
  data <- data.frame(p.value = pvalues, effectsize = effect_sizes, covariate1 = covariate, covariate2 = covariate2, is_alt = is_alt)
  return(data)
}

# Simulate data
num_tests <- 10000
df <- simulate_informative_data(num_tests)

# Linear regression model
model <- lm(effectsize ~ covariate1 + covariate2, data = df)

df$composite_covariate <- predict(model, newdata = df)

# Function to apply both IHW and BH on a subset of data and calculate discoveries
apply_methods <- function(subset_df) {
  pvalues = subset_df$p.value
  cov = subset_df$composite_covariate
  
  # Apply IHW
  ihw_result <- ihw(pvalues, cov, alpha = 0.05)
  discoveries_ihw <- sum(ihw_result@df$adj_pvalue <= 0.05)
  
  # Apply BH
  pvalues_adj_bh <- p.adjust(pvalues, method = "BH")
  discoveries_bh <- sum(pvalues_adj_bh <= 0.05)
  
  # Return the discoveries for both methods
  return(list(discoveries_ihw = discoveries_ihw, discoveries_bh = discoveries_bh))
}

# Function to create subsets and apply both methods
analyze_subsets <- function(df, subset_sizes) {
  results <- list()
  
  for (size in subset_sizes) {
    # Sample `size` rows from the dataframe
    subset_df <- df[sample(nrow(df), size), ]
    
    # Apply IHW and BH to the subset
    methods_result <- apply_methods(subset_df)
    
    # Store the result
    results[[paste("Size", size)]] <- c(methods_result, SubsetSize = size)
  }
  
  return(results)
}

subset_sizes <- seq(from = 3000, to = 10000, length.out = 8)

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
  labs(title = "Number of Discoveries vs. Subset Size with Regression",
       x = "Subset Size",
       y = "Number of Discoveries",
       color = "Method") +
  theme_minimal()
