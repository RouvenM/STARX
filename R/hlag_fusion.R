# hlag_fusion.R
# HLag + Kettenfusion fuer sparse VAR-Schaetzung.
# Erklaerungen, Herkunft der Funktionen und Referenzen siehe README.md.

prox_hlag_vec <- function(x, lambda) {
  p <- length(x)
  r <- x
  for (l in p:1) {
    idx <- l:p
    nr <- sqrt(sum(r[idx]^2))
    if (nr > 0) {
      shrink <- max(0, 1 - lambda / nr)
      r[idx] <- shrink * r[idx]
    }
  }
  r
}

prox_hlag_row <- function(Phi_i, lambda1) {
  apply(Phi_i, 2, prox_hlag_vec, lambda = lambda1)
}

build_diff_operator <- function(p) {
  D <- matrix(0, p - 1, p)
  for (l in 1:(p - 1)) {
    D[l, l] <- -1
    D[l, l + 1] <- 1
  }
  D
}

prox_fusion_row <- function(V, lambda2, weights, D, rho = 1,
                             n_iter = 150, tol = 1e-6) {
  p <- nrow(V); k <- ncol(V)
  if (p == 1 || lambda2 == 0) return(V)
  M <- diag(p) + rho * crossprod(D)
  Minv <- solve(M)
  Z <- V
  W <- D %*% Z
  U <- matrix(0, p - 1, k)
  for (it in 1:n_iter) {
    Z_new <- Minv %*% (V + rho * crossprod(D, W - U))
    DZ <- D %*% Z_new
    target <- DZ + U
    W_new <- matrix(0, p - 1, k)
    for (l in 1:(p - 1)) {
      v <- target[l, ]
      nv <- sqrt(sum(v^2))
      thresh <- lambda2 * weights[l] / rho
      W_new[l, ] <- if (nv > thresh) (1 - thresh / nv) * v else rep(0, k)
    }
    U <- U + DZ - W_new
    if (max(abs(Z_new - Z)) < tol) { Z <- Z_new; W <- W_new; break }
    Z <- Z_new; W <- W_new
  }
  Z
}

ridge_prelim_row <- function(Zfull, y_i, p, k, ridge_lambda = 1) {
  Teff <- length(y_i)
  A <- crossprod(Zfull) / Teff + ridge_lambda * diag(ncol(Zfull))
  b <- crossprod(Zfull, y_i) / Teff
  sol <- solve(A, b)
  t(matrix(sol, nrow = k, ncol = p))
}

adaptive_fusion_weights <- function(Phi_prelim, alpha_decay = 1, weight_eps = 1e-2) {
  p <- nrow(Phi_prelim)
  sapply(1:(p - 1), function(l) {
    dnorm_l <- sqrt(sum((Phi_prelim[l, ] - Phi_prelim[l + 1, ])^2))
    exp(-alpha_decay) / (dnorm_l + weight_eps)
  })
}

prox_combined_row <- function(V, lambda1, lambda2, weights, D,
                               dykstra_iter = 30, fusion_rho = 1) {
  x <- V
  p_corr <- matrix(0, nrow(V), ncol(V))
  q_corr <- matrix(0, nrow(V), ncol(V))
  for (it in 1:dykstra_iter) {
    y <- prox_hlag_row(x + p_corr, lambda1)
    p_corr <- x + p_corr - y
    x_new <- prox_fusion_row(y + q_corr, lambda2, weights, D, rho = fusion_rho)
    q_corr <- y + q_corr - x_new
    if (max(abs(x_new - x)) < 1e-7) { x <- x_new; break }
    x <- x_new
  }
  x
}

fit_hlag_fusion_var <- function(Y, p, lambda1, lambda2, alpha_decay = 1,
                                 adaptive = FALSE, ridge_lambda = 1, weight_eps = 1e-2,
                                 max_iter = 500, tol = 1e-6, verbose = FALSE) {
  Tn <- nrow(Y); k <- ncol(Y)
  Teff <- Tn - p
  Zlags <- lapply(1:p, function(l) Y[(p - l + 1):(Tn - l), , drop = FALSE])
  Ymat <- Y[(p + 1):Tn, , drop = FALSE]
  Zfull <- do.call(cbind, Zlags)

  L <- max(eigen(crossprod(Zfull) / Teff, symmetric = TRUE, only.values = TRUE)$values)
  step <- 1 / L

  D <- build_diff_operator(p)
  weights_const <- rep(exp(-alpha_decay), p - 1)

  Phi <- array(0, dim = c(p, k, k))
  weights_used <- matrix(NA, p - 1, k)

  for (i in 1:k) {
    if (adaptive) {
      Phi_prelim <- ridge_prelim_row(Zfull, Ymat[, i], p, k, ridge_lambda)
      weights_i <- adaptive_fusion_weights(Phi_prelim, alpha_decay, weight_eps)
    } else {
      weights_i <- weights_const
    }
    weights_used[, i] <- weights_i

    Phi_i <- matrix(0, p, k)
    theta <- Phi_i
    t_old <- 1
    it_used <- max_iter
    for (it in 1:max_iter) {
      pred <- Reduce(`+`, lapply(1:p, function(l) Zlags[[l]] %*% theta[l, ]))
      resid <- Ymat[, i] - pred
      grad <- t(sapply(1:p, function(l) crossprod(Zlags[[l]], resid) / Teff))
      V <- theta + step * grad
      Phi_new <- prox_combined_row(V, step * lambda1, step * lambda2, weights_i, D)
      t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
      theta <- Phi_new + ((t_old - 1) / t_new) * (Phi_new - Phi_i)
      if (max(abs(Phi_new - Phi_i)) < tol) { Phi_i <- Phi_new; it_used <- it; break }
      Phi_i <- Phi_new; t_old <- t_new
    }
    Phi[, , i] <- Phi_i
    if (verbose) {
      cat("Gleichung", i, "konvergiert nach", it_used, "Iterationen",
          if (adaptive) paste0(" | Gewichte: ", paste(round(weights_i, 2), collapse = ", ")) else "",
          "\n")
    }
  }
  attr(Phi, "weights_used") <- weights_used
  Phi
}
