# Custom bigtime Functions

This repository contains modified R and C++ functions adapted from the `bigtime` package, specifically featuring custom implementation for the `LagFused` penalty.

## Prerequisites

Make sure you have the following R packages installed:

```R
install.packages(c("Rcpp", "RcppArmadillo", "devtools"))
```

## How to Load the Functions

Source all required functions directly into your R session via:

```R
# Load required libraries
library(Rcpp)
library(devtools)

# Load all required R and C++ functions 
source_url("https://raw.githubusercontent.com/RouvenM/STARX/main/R/load_all.R")
```

## Reproducible Examples (Testing LagFused Penalty)

Below is a self-contained test script verifying row fusion under the `LagFused` penalty for `VAR(2)` and `VAR(4)` models.

```R
# Fix test data (3 standardized variables, n = 100)
set.seed(123)
n <- 100
k <- 3
p_p2 <- 2
Y <- scale(matrix(rnorm(n * k), ncol = k))
colnames(Y) <- c("RegionA", "RegionB", "RegionC")

# Helper function to compute absolute differences between Lag 1 and Lag 2
calc_lag_diffs <- function(Phi_matrix, k) {
  Lag1 <- Phi_matrix[, 1:k]
  Lag2 <- Phi_matrix[, (k + 1):(2 * k)]
  abs(Lag2 - Lag1)
}

# --- PART 1: Test for VAR(2) (p = 2) ---

# Step 1: Weak Lambda (lambda = 0.05)
message("\n--- PART 1 (p = 2) - STEP 1: Weak Lambda (lambda = 0.05) ---")
fit1_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "none", VARlseq = 0.05)
diffs1_p2 <- calc_lag_diffs(fit1_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs1_p2, 5))

# Step 2: Medium Lambda (lambda = 0.5) & Comparison with Step 1
message("\n--- PART 1 (p = 2) - STEP 2: Medium Lambda (lambda = 0.5) ---")
fit2_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "none", VARlseq = 0.5)
diffs2_p2 <- calc_lag_diffs(fit2_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs2_p2, 5))

# Verification check
is_smaller_or_equal_p2 <- all(diffs2_p2 <= diffs1_p2 + 1e-8)
if (is_smaller_or_equal_p2) {
  message("SUCCESS (p=2): All lag differences decreased or stayed equal compared to Step 1!")
} else {
  warning("CHECK FAILED (p=2): Some lag differences unexpectedly increased.")
}

# Step 3: Strong Lambda (lambda = 6.0) & Check for partial fusions
message("\n--- PART 1 (p = 2) - STEP 3: Strong Lambda (lambda = 6.0) ---")
fit3_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "none", VARlseq = 6.0)
diffs3_p2 <- calc_lag_diffs(fit3_p2$Phi[,,1], k)

rows_fused_step3_p2 <- apply(diffs3_p2, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p2[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda (lambda = 100.0) & Global fused check
message("\n--- PART 1 (p = 2) - STEP 4: Extremely Strong Lambda (lambda = 100.0) ---")
fit4_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "none", VARlseq = 100.0)
diffs4_p2 <- calc_lag_diffs(fit4_p2$Phi[,,1], k)

rows_fused_step4_p2 <- apply(diffs4_p2, 1, function(row) all(row < 1e-5))
all_rows_fused_p2 <- all(rows_fused_step4_p2)

if (all_rows_fused_p2) {
  message("SUCCESS (p=2): All rows have been completely and exactly fused under an extremely strong penalty!")
} else {
  warning("CHECK FAILED (p=2): Even with lambda = 100, not all rows were fused.")
}

# Step 5: Cross-Validation & Visualization
message("\n--- PART 1 (p = 2) - STEP 5: Cross-Validation & Visualization ---")
fit_cv_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p2$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p2$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p2)


# --- PART 2: Generalized Cascade Test for VAR(4) (p = 4) ---

p_p4 <- 4

# Helper function to compute maximum absolute differences across all adjacent lags (l vs. l-1)
calc_lag_diffs_max <- function(Phi_matrix, k = NULL, p = NULL, lambda_idx = NULL) {
  if (length(dim(Phi_matrix)) == 3) {
    n_lambda <- dim(Phi_matrix)[3]
    if (is.null(lambda_idx)) lambda_idx <- n_lambda
    Phi_matrix <- Phi_matrix[, , lambda_idx, drop = TRUE]
  }
  
  if (is.null(k)) k <- nrow(Phi_matrix)
  if (is.null(p)) p <- ncol(Phi_matrix) / k
  
  if (p < 2) stop(sprintf("For lag differences, p must be >= 2 (currently p = %.1f).", p))
  if (ncol(Phi_matrix) %% k != 0) {
    stop(sprintf("Number of columns (%d) is not a multiple of k (%d).", ncol(Phi_matrix), k))
  }
  
  max_diffs <- matrix(0, nrow = k, ncol = k)
  
  for (l in 2:p) {
    cols_prev <- ((l - 2) * k + 1):((l - 1) * k)
    cols_curr <- ((l - 1) * k + 1):(l * k)
    
    Lag_prev <- Phi_matrix[, cols_prev, drop = FALSE]
    Lag_curr <- Phi_matrix[, cols_curr, drop = FALSE]
    
    diff_curr <- abs(Lag_curr - Lag_prev)
    max_diffs <- pmax(max_diffs, diff_curr)
  }
  
  colnames(max_diffs) <- paste0("Var_", 1:k)
  rownames(max_diffs) <- paste0("Var_", 1:k)
  return(max_diffs)
}

# Step 1: Weak Lambda (lambda = 0.05)
message("\n--- PART 2 (p = 4) - STEP 1: Weak Lambda (lambda = 0.05) ---")
fit1_p4 <- sparseVAR(Y = Y, p = p_p4, VARpen = "LagFused", selection = "none", VARlseq = 0.05)
diffs1_p4 <- calc_lag_diffs_max(fit1_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs1_p4, 5))

# Step 2: Medium Lambda (lambda = 0.5) & Comparison with Step 1
message("\n--- PART 2 (p = 4) - STEP 2: Medium Lambda (lambda = 0.5) ---")
fit2_p4 <- sparseVAR(Y = Y, p = p_p4, VARpen = "LagFused", selection = "none", VARlseq = 0.5)
diffs2_p4 <- calc_lag_diffs_max(fit2_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs2_p4, 5))

is_smaller_or_equal_p4 <- all(diffs2_p4 <= diffs1_p4 + 1e-8)
if (is_smaller_or_equal_p4) {
  message("SUCCESS (p=4): All maximum lag differences decreased or stayed equal compared to Step 1!")
} else {
  warning("CHECK FAILED (p=4): Some lag differences unexpectedly increased.")
}

# Step 3: Strong Lambda (lambda = 30000.0) & Check for partial fusions
message("\n--- PART 2 (p = 4) - STEP 3: Strong Lambda (lambda = 30000.0) ---")
fit3_p4 <- sparseVAR(Y = Y, p = p_p4, VARpen = "LagFused", selection = "none", VARlseq = 30000.0)
diffs3_p4 <- calc_lag_diffs_max(fit3_p4$Phi)

rows_fused_step3_p4 <- apply(diffs3_p4, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p4[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused across all lags.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused across all lags.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda (lambda = 100000.0) & Global fused check
message("\n--- PART 2 (p = 4) - STEP 4: Extremely Strong Lambda (lambda = 100000.0) ---")
fit4_p4 <- sparseVAR(Y = Y, p = p_p4, VARpen = "LagFused", selection = "none", VARlseq = 100000.0)
diffs4_p4 <- calc_lag_diffs_max(fit4_p4$Phi)

rows_fused_step4_p4 <- apply(diffs4_p4, 1, function(row) all(row < 1e-5))
all_rows_fused_p4 <- all(rows_fused_step4_p4)

if (all_rows_fused_p4) {
  message("SUCCESS (p=4): All rows have been completely and exactly fused across all 4 lags!")
} else {
  warning("CHECK FAILED (p=4): Even with lambda = 100000, not all rows were completely fused.")
}

# Step 5: Cross-Validation & Visualization
message("\n--- PART 2 (p = 4) - STEP 5: Cross-Validation & Visualization ---")
fit_cv_p4 <- sparseVAR(Y = Y, p = p_p4, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p4$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p4$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p4)
```




