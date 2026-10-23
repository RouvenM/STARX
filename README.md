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

## Reproducible Examples (Testing LagFused Penalty) on simulated white noise

Below is a self-contained test script verifying row fusion under the `LagFused` penalty for `VAR(2)` and `VAR(4)` models.

```R
# Fix test data (3 standardized variables, n = 2000)
set.seed(123)
n <- 2000
k <- 3
p <- 2
Y <- scale(matrix(rnorm(n * k), ncol = k))
colnames(Y) <- c("RegionA", "RegionB", "RegionC")

# Helper function to compute absolute differences between Lag 1 and Lag 2
calc_lag_diffs <- function(Phi_matrix, k) {
  Lag1 <- Phi_matrix[, 1:k]
  Lag2 <- Phi_matrix[, (k + 1):(2 * k)]
  abs(Lag2 - Lag1)
}

# --- PART 1: Test for VAR(2) (p = 2) ---

# --- Lambda values ---
lambda_weak    <- 0.05
lambda_medium  <- 0.5
lambda_strong  <- 8.0
lambda_extreme <- 100.0

# Step 1: Weak Lambda
message(sprintf("\n--- PART 1 (p = %d) - STEP 1: Weak Lambda (lambda = %g) ---", p, lambda_weak))
fit1_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_weak)
diffs1_p2 <- calc_lag_diffs(fit1_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs1_p2, 5))

# Step 2: Medium Lambda & Comparison with Step 1
message(sprintf("\n--- PART 1 (p = %d) - STEP 2: Medium Lambda (lambda = %g) ---", p, lambda_medium))
fit2_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_medium)
diffs2_p2 <- calc_lag_diffs(fit2_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs2_p2, 5))

# Verification check
is_smaller_or_equal_p2 <- all(diffs2_p2 <= diffs1_p2 + 1e-8)
if (is_smaller_or_equal_p2) {
  message(sprintf("SUCCESS (p=%d): All lag differences decreased or stayed equal compared to Step 1!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Some lag differences unexpectedly increased.", p))
}

# Step 3: Strong Lambda & Check for partial fusions
message(sprintf("\n--- PART 1 (p = %d) - STEP 3: Strong Lambda (lambda = %g) ---", p, lambda_strong))
fit3_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_strong)
diffs3_p2 <- calc_lag_diffs(fit3_p2$Phi[,,1], k)

rows_fused_step3_p2 <- apply(diffs3_p2, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p2[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda & Global fused check
message(sprintf("\n--- PART 1 (p = %d) - STEP 4: Extremely Strong Lambda (lambda = %g) ---", p, lambda_extreme))
fit4_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_extreme)
diffs4_p2 <- calc_lag_diffs(fit4_p2$Phi[,,1], k)

rows_fused_step4_p2 <- apply(diffs4_p2, 1, function(row) all(row < 1e-5))
all_rows_fused_p2 <- all(rows_fused_step4_p2)

if (all_rows_fused_p2) {
  message(sprintf("SUCCESS (p=%d): All rows have been completely and exactly fused under an extremely strong penalty!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Even with lambda = %g, not all rows were fused.", p, lambda_extreme))
}

# Step 5: Cross-Validation & Visualization
message(sprintf("\n--- PART 1 (p = %d) - STEP 5: Cross-Validation & Visualization ---", p))
fit_cv_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p2$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p2$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p2)


# --- PART 2: Generalized Cascade Test for VAR(4) (p = 4) ---

p <- 4

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

# --- PART 2: Test for VAR(4) (p = 4) ---

# --- Lambda values ---
lambda_weak    <- 0.05
lambda_medium  <- 0.5
lambda_strong  <- 220471.0
lambda_extreme <- 300000.0

# Step 1: Weak Lambda
message(sprintf("\n--- PART 2 (p = %d) - STEP 1: Weak Lambda (lambda = %g) ---", p, lambda_weak))
fit1_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_weak)
diffs1_p4 <- calc_lag_diffs_max(fit1_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs1_p4, 5))

# Step 2: Medium Lambda & Comparison with Step 1
message(sprintf("\n--- PART 2 (p = %d) - STEP 2: Medium Lambda (lambda = %g) ---", p, lambda_medium))
fit2_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_medium)
diffs2_p4 <- calc_lag_diffs_max(fit2_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs2_p4, 5))

is_smaller_or_equal_p4 <- all(diffs2_p4 <= diffs1_p4 + 1e-8)
if (is_smaller_or_equal_p4) {
  message(sprintf("SUCCESS (p=%d): All maximum lag differences decreased or stayed equal compared to Step 1!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Some lag differences unexpectedly increased.", p))
}

# Step 3: Strong Lambda & Check for partial fusions
message(sprintf("\n--- PART 2 (p = %d) - STEP 3: Strong Lambda (lambda = %g) ---", p, lambda_strong))
fit3_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_strong)
diffs3_p4 <- calc_lag_diffs_max(fit3_p4$Phi)

rows_fused_step3_p4 <- apply(diffs3_p4, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p4[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused across all lags.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused across all lags.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda & Global fused check
message(sprintf("\n--- PART 2 (p = %d) - STEP 4: Extremely Strong Lambda (lambda = %g) ---", p, lambda_extreme))
fit4_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_extreme)
diffs4_p4 <- calc_lag_diffs_max(fit4_p4$Phi)

rows_fused_step4_p4 <- apply(diffs4_p4, 1, function(row) all(row < 1e-5))
all_rows_fused_p4 <- all(rows_fused_step4_p4)

if (all_rows_fused_p4) {
  message(sprintf("SUCCESS (p=%d): All rows have been completely and exactly fused across all %d lags!", p, p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Even with lambda = %g, not all rows were completely fused.", p, lambda_extreme))
}

# Step 5: Cross-Validation & Visualization
message(sprintf("\n--- PART 2 (p = %d) - STEP 5: Cross-Validation & Visualization ---", p))
fit_cv_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p4$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p4$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p4)
```

## Reproducible Examples (Testing LagFused Penalty) on simulated VAR processes

Below is a self-contained test script verifying row fusion under the `LagFused` penalty for `VAR(2)` and `VAR(4)` models.

```R
# --- PART 1: Test for VAR(2) (p = 2) ---

k <- 3
p <- 2

# Define lag matrices (Row 1 is identical across lags / fused)
Phi1 <- matrix(c(
  0.4,  0.2,  0.0,  # Region A (fused)
  0.3, -0.2,  0.1,  # Region B
  0.0,  0.1,  0.4   # Region C
), nrow = k, byrow = TRUE)

Phi2 <- matrix(c(
  0.4,  0.2,  0.0,  # Region A (fused with Phi1[1, ])
 -0.2,  0.0,  0.1,  # Region B
  0.0,  0.0, -0.2   # Region C
), nrow = k, byrow = TRUE)

# Construct coefficient matrix in companion form (kp x kp)
top_block <- cbind(Phi1, Phi2)
bottom_block <- cbind(diag(k), matrix(0, nrow = k, ncol = k))
coef_mat_companion <- rbind(top_block, bottom_block)

# Verify stationarity (modulus of all eigenvalues must be < 1)
max_ev <- max(abs(eigen(coef_mat_companion)$values))
if (max_ev >= 1) {
  stop(paste("Process is not stationary! Max eigenvalue:", round(max_ev, 3)))
}

# Simulate VAR(2) process
sim_res <- simVAR(
  periods = 2000,
  k = k,
  p = p,
  coef_mat = coef_mat_companion, 
  seed = 123
)

# Extract and scale simulated time series
Y <- scale(sim_res$Y)
colnames(Y) <- c("RegionA", "RegionB", "RegionC")

# --- Lambda values ---
lambda_weak    <- 0.05
lambda_medium  <- 0.5
lambda_strong  <- 100.0
lambda_extreme <- 1000.0

# Step 1: Weak Lambda
message(sprintf("\n--- PART 1 (p = %d) - STEP 1: Weak Lambda (lambda = %g) ---", p, lambda_weak))
fit1_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_weak)
diffs1_p2 <- calc_lag_diffs(fit1_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs1_p2, 5))

# Step 2: Medium Lambda & Comparison with Step 1
message(sprintf("\n--- PART 1 (p = %d) - STEP 2: Medium Lambda (lambda = %g) ---", p, lambda_medium))
fit2_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_medium)
diffs2_p2 <- calc_lag_diffs(fit2_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs2_p2, 5))

# Verification check
is_smaller_or_equal_p2 <- all(diffs2_p2 <= diffs1_p2 + 1e-8)
if (is_smaller_or_equal_p2) {
  message(sprintf("SUCCESS (p=%d): All lag differences decreased or stayed equal compared to Step 1!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Some lag differences unexpectedly increased.", p))
}

# Step 3: Strong Lambda & Check for partial fusions
message(sprintf("\n--- PART 1 (p = %d) - STEP 3: Strong Lambda (lambda = %g) ---", p, lambda_strong))
fit3_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_strong)
diffs3_p2 <- calc_lag_diffs(fit3_p2$Phi[,,1], k)

rows_fused_step3_p2 <- apply(diffs3_p2, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p2[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda & Global fused check
message(sprintf("\n--- PART 1 (p = %d) - STEP 4: Extremely Strong Lambda (lambda = %g) ---", p, lambda_extreme))
fit4_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_extreme)
diffs4_p2 <- calc_lag_diffs(fit4_p2$Phi[,,1], k)

rows_fused_step4_p2 <- apply(diffs4_p2, 1, function(row) all(row < 1e-5))
all_rows_fused_p2 <- all(rows_fused_step4_p2)

if (all_rows_fused_p2) {
  message(sprintf("SUCCESS (p=%d): All rows have been completely and exactly fused under an extremely strong penalty!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Even with lambda = %g, not all rows were fused.", p, lambda_extreme))
}

# Step 5: Cross-Validation & Visualization
message(sprintf("\n--- PART 1 (p = %d) - STEP 5: Cross-Validation & Visualization ---", p))
fit_cv_p2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p2$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p2$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p2)

# --- PART 2: Generalized Cascade Test for VAR(4) (p = 4) ---

k <- 3
p <- 4

# Define row 1 values (identical across all 4 lag matrices)
row1_fused <- c(0.15, 0.05, 0.00)  # Region A: completely fused

# 1. Define 4 lag matrices (3x3 each)
Phi1 <- matrix(c(
  row1_fused,           # Region A (fused)
  0.20, -0.10,  0.10,   # Region B
  0.00,  0.10,  0.20    # Region C
), nrow = k, byrow = TRUE)

Phi2 <- matrix(c(
  row1_fused,           # Region A (fused)
 -0.10,  0.10,  0.00,   # Region B
  0.10,  0.00, -0.10    # Region C
), nrow = k, byrow = TRUE)

Phi3 <- matrix(c(
  row1_fused,           # Region A (fused)
  0.00, -0.15,  0.10,   # Region B
 -0.10,  0.10,  0.00    # Region C
), nrow = k, byrow = TRUE)

Phi4 <- matrix(c(
  row1_fused,           # Region A (fused)
  0.10,  0.00, -0.10,   # Region B
  0.00,  0.00,  0.10    # Region C
), nrow = k, byrow = TRUE)

# 2. Construct general companion matrix for p = 4 (12x12)
top_block <- cbind(Phi1, Phi2, Phi3, Phi4)
bottom_block <- cbind(diag((p - 1) * k), matrix(0, nrow = (p - 1) * k, ncol = k))
coef_mat_companion <- rbind(top_block, bottom_block)

# 3. Verify stationarity
max_ev <- max(abs(eigen(coef_mat_companion)$values))
if (max_ev >= 1) {
  stop(paste("Process is not stationary! Max eigenvalue:", round(max_ev, 3)))
}

# 4. Simulate VAR(4) process
sim_res <- simVAR(
  k = k,
  p = p,
  periods = 2000,
  coef_mat = coef_mat_companion,
  seed = 123
)

# 5. Extract and scale simulated time series
Y <- scale(sim_res$Y)
colnames(Y) <- c("RegionA", "RegionB", "RegionC")

# --- PART 2: Test for VAR(4) (p = 4) ---

# --- Lambda values ---
lambda_weak    <- 0.05
lambda_medium  <- 0.5
lambda_strong  <- 500000.0
lambda_extreme <- 1000000.0

# Step 1: Weak Lambda
message(sprintf("\n--- PART 2 (p = %d) - STEP 1: Weak Lambda (lambda = %g) ---", p, lambda_weak))
fit1_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_weak)
diffs1_p4 <- calc_lag_diffs_max(fit1_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs1_p4, 5))

# Step 2: Medium Lambda & Comparison with Step 1
message(sprintf("\n--- PART 2 (p = %d) - STEP 2: Medium Lambda (lambda = %g) ---", p, lambda_medium))
fit2_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_medium)
diffs2_p4 <- calc_lag_diffs_max(fit2_p4$Phi)

cat("Maximum absolute differences max|Phi_ij^(l) - Phi_ij^(l-1)|:\n")
print(round(diffs2_p4, 5))

is_smaller_or_equal_p4 <- all(diffs2_p4 <= diffs1_p4 + 1e-8)
if (is_smaller_or_equal_p4) {
  message(sprintf("SUCCESS (p=%d): All maximum lag differences decreased or stayed equal compared to Step 1!", p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Some lag differences unexpectedly increased.", p))
}

# Step 3: Strong Lambda & Check for partial fusions
message(sprintf("\n--- PART 2 (p = %d) - STEP 3: Strong Lambda (lambda = %g) ---", p, lambda_strong))
fit3_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_strong)
diffs3_p4 <- calc_lag_diffs_max(fit3_p4$Phi)

rows_fused_step3_p4 <- apply(diffs3_p4, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3_p4[i]) {
    cat(sprintf("  -> Row %d (%s) IS exactly fused across all lags.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Row %d (%s) is NOT YET completely fused across all lags.\n", i, colnames(Y)[i]))
  }
}

# Step 4: Extremely Strong Lambda & Global fused check
message(sprintf("\n--- PART 2 (p = %d) - STEP 4: Extremely Strong Lambda (lambda = %g) ---", p, lambda_extreme))
fit4_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = lambda_extreme)
diffs4_p4 <- calc_lag_diffs_max(fit4_p4$Phi)

rows_fused_step4_p4 <- apply(diffs4_p4, 1, function(row) all(row < 1e-5))
all_rows_fused_p4 <- all(rows_fused_step4_p4)

if (all_rows_fused_p4) {
  message(sprintf("SUCCESS (p=%d): All rows have been completely and exactly fused across all %d lags!", p, p))
} else {
  warning(sprintf("CHECK FAILED (p=%d): Even with lambda = %g, not all rows were completely fused.", p, lambda_extreme))
}

# Step 5: Cross-Validation & Visualization
message(sprintf("\n--- PART 2 (p = %d) - STEP 5: Cross-Validation & Visualization ---", p))
fit_cv_p4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimal Lambda (Min MSFE)  : %.4f\n", fit_cv_p4$lambda_opt))
cat(sprintf("  -> Optimal Lambda (1-SE Rule) : %.4f\n", fit_cv_p4$lambda_SEopt))

# Visualization: Cross-Validation error curve
plot_cv(fit_cv_p4)
```





