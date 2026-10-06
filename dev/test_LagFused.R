# ==============================================================================
# IMB: Stufenweiser Test der Time-Lag Fusion über verschiedene Lambda-Ränge
# ==============================================================================

devtools::document()
devtools::load_all()

# Fixierte Testdaten (VAR(p) mit 3 Variablen, standardisiert)
set.seed(123)
n <- 100
k <- 3
p <- 4
# p <- 2
Y <- scale(matrix(rnorm(n * k), ncol = k))
colnames(Y) <- c("RegionA", "RegionB", "RegionC")

# Hilfsfunktion zur Berechung der absoluten Abstände zwischen Lag 1 und Lag 2
calc_lag_diffs <- function(Phi_matrix, k) {
  Lag1 <- Phi_matrix[, 1:k]
  Lag2 <- Phi_matrix[, (k + 1):(2 * k)]
  abs(Lag2 - Lag1)
}


# ------------------------------------------------------------------------------
# Schritt 1: Schwaches Lambda (lambda = 0.05)
# ------------------------------------------------------------------------------
message("\n=== SCHRITT 1: Schwaches Lambda (lambda = 0.05) ===")
fit1 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = 0.05)
diffs1 <- calc_lag_diffs(fit1$Phi[,,1], k)

cat("Absolute Abstände |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs1, 5))


# ------------------------------------------------------------------------------
# Schritt 2: Mittleres Lambda (lambda = 0.5) & Vergleich zu Schritt 1
# ------------------------------------------------------------------------------
message("\n=== SCHRITT 2: Mittleres Lambda (lambda = 0.5) ===")
fit2 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = 0.5)
diffs2 <- calc_lag_diffs(fit2$Phi[,,1], k)

cat("Absolute Abstände |Phi_ij^(2) - Phi_ij^(1)|:\n")
print(round(diffs2, 5))

is_smaller_or_equal <- all(diffs2 <= diffs1 + 1e-8)
if (is_smaller_or_equal) {
  message("SUCCESS: Alle Abstände sind im Vergleich zu Schritt 1 kleiner oder gleich geworden!")
} else {
  warning("CHECK FAILED: Einige Abstände sind unerwartet gewachsen.")
}


# ------------------------------------------------------------------------------
# Schritt 3: Starkes Lambda (lambda = 6.0) & Prüfung partieller Fusionen
# ------------------------------------------------------------------------------
message("\n=== SCHRITT 3: Starkes Lambda (lambda = 6.0) ===")
fit3 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = 6.0)
diffs3 <- calc_lag_diffs(fit3$Phi[,,1], k)

rows_fused_step3 <- apply(diffs3, 1, function(row) all(row < 1e-5))
for (i in 1:k) {
  if (rows_fused_step3[i]) {
    cat(sprintf("  -> Zeile %d (%s) IST exakt fusioniert.\n", i, colnames(Y)[i]))
  } else {
    cat(sprintf("  -> Zeile %d (%s) ist noch NICHT vollständig fusioniert.\n", i, colnames(Y)[i]))
  }
}


# ------------------------------------------------------------------------------
# Schritt 4: Extrem starkes Lambda (lambda = 100.0) & Globale Fused-Prüfung
# ------------------------------------------------------------------------------
message("\n=== SCHRITT 4: Extrem starkes Lambda (lambda = 100.0) ===")
fit4 <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "none", VARlseq = 100.0)
diffs4 <- calc_lag_diffs(fit4$Phi[,,1], k)

rows_fused_step4 <- apply(diffs4, 1, function(row) all(row < 1e-5))
all_rows_fused <- all(rows_fused_step4)

if (all_rows_fused) {
  message("SUCCESS: Alle Zeilen wurden bei extrem starker Strafe vollständig und exakt fusioniert!")
} else {
  warning("CHECK FAILED: Selbst bei lambda = 100 wurden nicht alle Zeilen fusioniert.")
}


# ------------------------------------------------------------------------------
# Schritt 5: Cross-Validation & Visualisierung
# ------------------------------------------------------------------------------
message("\n=== SCHRITT 5: Cross-Validation & Visualisierung ===")
fit_cv <- sparseVAR(Y = Y, p = p, VARpen = "LagFused", selection = "cv")

cat(sprintf("  -> Optimales Lambda (Min MSFE)  : %.4f\n", fit_cv$lambda_opt))
cat(sprintf("  -> Optimales Lambda (1-SE Rule) : %.4f\n", fit_cv$lambda_SEopt))

# Visualisierung 1: Cross-Validation Fehlerkurve
plot_cv(fit_cv)

# Visualisierung 2: Lag-Matrix Heatmap
lagmatrix(fit = fit_cv, returnplot = TRUE)






calc_lag_diffs_max <- function(Phi_matrix, k = NULL, p = NULL, lambda_idx = NULL) {
  # 1. 3D-Array Behandlung (falls Lambda-Gitter übergeben wird)
  if (length(dim(Phi_matrix)) == 3) {
    n_lambda <- dim(Phi_matrix)[3]
    if (is.null(lambda_idx)) lambda_idx <- n_lambda
    Phi_matrix <- Phi_matrix[, , lambda_idx, drop = TRUE]
  }
  
  # 2. Automatische Bestimmung von k und p (falls nicht angegeben)
  if (is.null(k)) k <- nrow(Phi_matrix)
  if (is.null(p)) p <- ncol(Phi_matrix) / k
  
  if (p < 2) stop(sprintf("Für Lag-Differenzen muss p >= 2 sein (aktuell p = %.1f).", p))
  if (ncol(Phi_matrix) %% k != 0) {
    stop(sprintf("Spaltenanzahl (%d) ist kein Vielfaches von k (%d).", ncol(Phi_matrix), k))
  }
  
  # Matrix für maximale Abstände initialisieren (k x k)
  max_diffs <- matrix(0, nrow = k, ncol = k)
  
  # Schleife über alle benachbarten Lag-Blöcke (l vs. l-1)
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

# Berechnet automatisch k = 3 und p = 4 aus der Matrixform!
diffs_cv <- calc_lag_diffs_max(fit_cv$Phihat)
print(round(diffs_cv, 5))

print(round(calc_lag_diffs_max(fit4$Phihat), 5))

fit_tight <- sparseVAR(
  Y = Y, p = p, VARpen = "LagFused", 
  selection = "none", VARlseq = 100000, 
  eps = 1e-12
)

round(calc_lag_diffs_max(fit_tight$Phihat), 5)

fit1$Phihat
fit2$Phihat
fit3$Phihat
fit4$Phihat

round(calc_lag_diffs_max(fit1$Phihat), 5)
round(calc_lag_diffs_max(fit2$Phihat), 5)
round(calc_lag_diffs_max(fit3$Phihat), 5)
round(calc_lag_diffs_max(fit4$Phihat), 5)






