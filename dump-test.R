# Example Script: Using NNEBP with Dump Files (Minimal Version)

library(povmap)
library(ggplot2)
library(dplyr)

# Source the updated NNEBP implementation
source("nnebp_v6.R")

# Load example data
data("eusilcA_pop")
data("eusilcA_smp")

# Define formula - use a simple model for quick execution
formula <- eqIncome ~ gender + eqsize + cash + self_empl

# Define output file paths
ydump_file <- "household_level_simulations.csv"
pdump_file <- "area_level_simulations.csv"

# Run NNEBP with dump files
nnebp_result <- nnebp(
  fixed = formula,
  pop_data = eusilcA_pop,
  pop_domains = "district",
  smp_data = eusilcA_smp,
  smp_domains = "district",
  L = 5,  # Small number for demonstration
  MSE = TRUE,
  B = 3,  # Small number for demonstration
  mse_method = "variance",
  na.rm = TRUE,
  Ydump = ydump_file,
  Pdump = pdump_file
)

# Read the area-level simulation data
area_data <- read.csv(pdump_file)

# Print column names to verify
cat("Column names in Pdump file:", paste(colnames(area_data), collapse = ", "), "\n")

# Verify column names actually match what we expect
expected_cols <- c("L", "Domain", "Area_Effect", "Mean", "Head_Count", 
                   "Poverty_Gap", "Gini", "Median", "Quantile_10", 
                   "Quantile_25", "Quantile_75", "Quantile_90", "Quintile_Share")
missing_cols <- setdiff(expected_cols, colnames(area_data))
if (length(missing_cols) > 0) {
  cat("Missing expected columns:", paste(missing_cols, collapse = ", "), "\n")
  
  # Try reading with check.names = FALSE
  area_data <- read.csv(pdump_file, check.names = FALSE)
  cat("Column names with check.names=FALSE:", paste(colnames(area_data), collapse = ", "), "\n")
}

# If column names are X1, X2, etc. (no headers were recognized), manually rename them
if (all(grepl("^X[0-9]+$", colnames(area_data)[1:3]))) {
  cat("CSV seems to have been read without proper headers, renaming columns manually\n")
  colnames(area_data) <- expected_cols
}

# Now try plotting with correct column names
plot1 <- ggplot(area_data, aes(x = Domain, y = Area_Effect)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of Area Effects by Domain",
       y = "Area Effect")
print(plot1)

# Plot relationship between Area_Effect and Mean
plot2 <- ggplot(area_data, aes(x = Area_Effect, y = Mean)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  labs(title = "Relationship between Area Effect and Mean Estimate",
       x = "Area Effect", 
       y = "Mean Estimate")
print(plot2)

# Clean up (optional)
# file.remove(ydump_file)
# file.remove(pdump_file)