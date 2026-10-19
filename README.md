# Custom bigtime Functions

This repository contains modified R and C++ functions adapted from the `bigtime` package, specifically featuring custom implementation for the `LagFused` penalty.

## Prerequisites

Make sure you have the following R packages installed before loading the functions:

```R
install.packages(c("bigtime", "Rcpp", "devtools"))
```

## How to Load the Functions

You can source the required R and C++ functions directly into your R session from this GitHub repository without installing a full R package:

```R
# 1. Load required libraries
#library(bigtime)
library(Rcpp)
library(devtools)

# 2. Source R helper functions from bigtime
source_url("[https://raw.githubusercontent.com/RouvenM/STARX/main/R/diagnostics.R](https://raw.githubusercontent.com/RouvenM/STARX/main/R/diagnostics.R)")

# 3. Compile and load C++ routines directly from GitHub
sourceCpp("[https://raw.githubusercontent.com/RouvenM/STARX/main/src/hvar.cpp](https://raw.githubusercontent.com/RouvenM/STARX/main/src/hvar.cpp)")

# 4. Load modified main R functions
source_url("[https://raw.githubusercontent.com/RouvenM/STARX/main/R/hvar.R](https://raw.githubusercontent.com/RouvenM/STARX/main/R/hvar.R)")
```

## Reproducible Example (Testing LagFused Penalty)

Below is a self-contained test script verifying that higher values of `lambda` lead to increased row fusion between coefficient matrices $\Phi^{(1)}$ and $\Phi^{(2)}$ in a $\text{VAR}(2)$ model.

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

# Step 1: Weak Lambda (lambda = 0.05)
message("\n=== PART 1 (p = 2) - STEP 1: Weak Lambda (lambda = 0.05) ===")
fit1_p2 <- sparseVAR(Y = Y, p = p_p2, VARpen = "LagFused", selection = "none", VARlseq = 0.05)
diffs1_p2 <- calc_lag_diffs(fit1_p2$Phi[,,1], k)

cat("Absolute differences |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs1_p2, 5))

# Step 2: Medium Lambda (lambda = 0.5) & Comparison with Step 1
message("\n=== PART 1 (p = 2) - STEP 2: Medium Lambda (lambda = 0.5) ===")
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
```