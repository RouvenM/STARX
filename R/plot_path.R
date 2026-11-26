#' Plot Frobenius Norm Paths of Lag Matrices Across Lambda
#'
#' Re-fits the sparse VAR model over the saved lambda grid and plots the Frobenius
#' norm of each lag matrix \eqn{\Phi^{(\ell)}} as a function of \eqn{\lambda}.
#'
#' @param fit A fitted object returned by \code{\link{sparseVAR}} (or \code{cv.sparseVAR})
#'   with class \code{"bigtime.VAR"}.
#' @param Y Optional numeric matrix of time series data. If \code{NULL}, uses \code{fit$Y}.
#' @param VARpen Optional character string specifying the VAR penalty. If \code{NULL},
#'   uses \code{fit$VARpen}.
#' @param log_lambda Logical; if \code{TRUE}, the x-axis displays \eqn{\log(\lambda)}. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed down to \code{\link{sparseVAR}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing the Frobenius norm paths across the lambda grid.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs theme_minimal
plot_frobenius_path <- function(fit, Y = NULL, VARpen = NULL, log_lambda = TRUE, ...) {
  
  if (is.null(Y)) {
    Y <- fit$Y
    if (is.null(Y)) {
      stop("Data matrix Y was not found in 'fit' object. Please provide Y explicitly.")
    }
  }
  
  if (is.null(VARpen)) {
    VARpen <- fit$VARpen
    if (is.null(VARpen)) {
      stop("VAR penalty was not found in 'fit' object. Please provide 'VARpen' explicitly.")
    }
  }
  
  lambda_grid <- fit$lambdas
  if (is.null(lambda_grid)) {
    stop("No lambda grid found in the provided fit object.")
  }
  
  p <- fit$p
  k <- fit$k
  
  # Refit across the lambda grid and compute the Frobenius norm of each lag matrix
  norm_list <- lapply(lambda_grid, function(l) {
    fit_l <- sparseVAR(
      Y = Y,
      p = p,
      VARpen = VARpen,
      VARlseq = l,
      selection = "none",
      ...
    )
    
    phi_mat <- fit_l$Phihat
    if (length(dim(phi_mat)) == 3) {
      phi_mat <- matrix(phi_mat[, , 1], nrow = k, ncol = k * p)
    }
    
    lag_norms <- numeric(p)
    for (lag in seq_len(p)) {
      if (is.list(phi_mat)) {
        phi_lag <- phi_mat[[lag]]
      } else if (is.matrix(phi_mat)) {
        col_indices <- ((lag - 1) * k + 1):(lag * k)
        phi_lag <- phi_mat[, col_indices, drop = FALSE]
      } else {
        stop("Unexpected structure in fit_l$Phihat.")
      }
      
      lag_norms[lag] <- norm(phi_lag, type = "F")
    }
    
    data.frame(
      lambda = l,
      lag = factor(paste0("Lag ", seq_len(p)), levels = paste0("Lag ", seq_len(p))),
      frobenius_norm = lag_norms
    )
  })
  
  plot_data <- do.call(rbind, norm_list)
  
  # Extract CV-selected lambdas if available
  opt_val <- if (!is.null(fit$lambda_opt) && length(fit$lambda_opt) == 1 && !is.na(fit$lambda_opt)) {
    fit$lambda_opt
  } else {
    NULL
  }
  se_opt_val <- if (!is.null(fit$lambda_SEopt) && length(fit$lambda_SEopt) == 1 && !is.na(fit$lambda_SEopt)) {
    fit$lambda_SEopt
  } else {
    NULL
  }
  
  if (log_lambda) {
    plot_data$x_val <- log(plot_data$lambda)
    x_label <- expression(log(lambda))
    if (!is.null(opt_val)) opt_val <- log(opt_val)
    if (!is.null(se_opt_val)) se_opt_val <- log(se_opt_val)
  } else {
    plot_data$x_val <- plot_data$lambda
    x_label <- expression(lambda)
  }
  
  p_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = frobenius_norm, color = lag)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = x_label,
      y = expression(paste("||", Phi^{(l)}, "||"[F])),
      color = "Lag Matrix",
      title = "Frobenius Norm Paths across Penalty Grid"
    ) +
    ggplot2::theme_minimal()
  
  # Add vertical reference lines for optimal (dashed) and 1-SE (dotted) lambda
  if (!is.null(opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = opt_val, linetype = "dashed", color = "black", linewidth = 1.1)
  }
  if (!is.null(se_opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = se_opt_val, linetype = "dotted", color = "black", linewidth = 1.1)
  }
  
  return(p_obj)
}


#' Plot Row-wise Frobenius Norm Paths of Lag Matrices Across Lambda
#'
#' Re-fits the sparse VAR model over the saved lambda grid and plots the row-wise
#' Frobenius (L2) norms for each row of each lag matrix \eqn{\Phi^{(\ell)}} as a function of \eqn{\lambda}.
#'
#' @param fit A fitted object returned by \code{\link{sparseVAR}} (or \code{cv.sparseVAR})
#'   with class \code{"bigtime.VAR"}.
#' @param Y Optional numeric matrix of time series data. If \code{NULL}, uses \code{fit$Y}.
#' @param VARpen Optional character string specifying the VAR penalty. If \code{NULL},
#'   uses \code{fit$VARpen}.
#' @param log_lambda Logical; if \code{TRUE}, the x-axis displays \eqn{\log(\lambda)}. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed down to \code{\link{sparseVAR}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing row-wise norm paths across the lambda grid, faceted by target variable.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs theme_minimal facet_wrap
plot_row_frobenius_path <- function(fit, Y = NULL, VARpen = NULL, log_lambda = TRUE, ...) {
  
  if (is.null(Y)) {
    Y <- fit$Y
    if (is.null(Y)) {
      stop("Data matrix Y was not found in 'fit' object. Please provide Y explicitly.")
    }
  }
  
  if (is.null(VARpen)) {
    VARpen <- fit$VARpen
    if (is.null(VARpen)) {
      stop("VAR penalty was not found in 'fit' object. Please provide 'VARpen' explicitly.")
    }
  }
  
  lambda_grid <- fit$lambdas
  if (is.null(lambda_grid)) {
    stop("No lambda grid found in the provided fit object.")
  }
  
  p <- fit$p
  k <- fit$k
  series_names <- if (!is.null(fit$series_names)) fit$series_names else paste0("Row ", seq_len(k))
  
  # Refit across the lambda grid and compute row-wise L2 norms for each lag matrix
  norm_list <- lapply(lambda_grid, function(l) {
    fit_l <- sparseVAR(
      Y = Y,
      p = p,
      VARpen = VARpen,
      VARlseq = l,
      selection = "none",
      ...
    )
    
    phi_mat <- fit_l$Phihat
    if (length(dim(phi_mat)) == 3) {
      phi_mat <- matrix(phi_mat[, , 1], nrow = k, ncol = k * p)
    }
    
    grid_rows <- list()
    for (lag in seq_len(p)) {
      if (is.list(phi_mat)) {
        phi_lag <- phi_mat[[lag]]
      } else if (is.matrix(phi_mat)) {
        col_indices <- ((lag - 1) * k + 1):(lag * k)
        phi_lag <- phi_mat[, col_indices, drop = FALSE]
      } else {
        stop("Unexpected structure in fit_l$Phihat.")
      }
      
      row_norms <- sqrt(rowSums(phi_lag^2))
      
      grid_rows[[lag]] <- data.frame(
        lambda = l,
        lag = factor(paste0("Lag ", lag), levels = paste0("Lag ", seq_len(p))),
        row_variable = factor(series_names, levels = series_names),
        row_norm = row_norms
      )
    }
    
    do.call(rbind, grid_rows)
  })
  
  plot_data <- do.call(rbind, norm_list)
  
  # Extract CV-selected lambdas if available
  opt_val <- if (!is.null(fit$lambda_opt) && length(fit$lambda_opt) == 1 && !is.na(fit$lambda_opt)) {
    fit$lambda_opt
  } else {
    NULL
  }
  se_opt_val <- if (!is.null(fit$lambda_SEopt) && length(fit$lambda_SEopt) == 1 && !is.na(fit$lambda_SEopt)) {
    fit$lambda_SEopt
  } else {
    NULL
  }
  
  if (log_lambda) {
    plot_data$x_val <- log(plot_data$lambda)
    x_label <- expression(log(lambda))
    if (!is.null(opt_val)) opt_val <- log(opt_val)
    if (!is.null(se_opt_val)) se_opt_val <- log(se_opt_val)
  } else {
    plot_data$x_val <- plot_data$lambda
    x_label <- expression(lambda)
  }
  
  p_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = row_norm, color = lag)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~ row_variable) +
    ggplot2::labs(
      x = x_label,
      y = expression(paste("||", Phi[i, "."], "||"[2]^{(l)})),
      color = "Lag",
      title = "Row-wise Frobenius Norm Paths across Penalty Grid"
    ) +
    ggplot2::theme_minimal()
  
  # Add vertical reference lines for optimal (dashed) and 1-SE (dotted) lambda
  if (!is.null(opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = opt_val, linetype = "dashed", color = "black", linewidth = 1.1)
  }
  if (!is.null(se_opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = se_opt_val, linetype = "dotted", color = "black", linewidth = 1.1)
  }
  
  return(p_obj)
}


#' Plot Individual Coefficient Paths Across Lambda
#'
#' Re-fits the sparse VAR model over the saved lambda grid and plots the raw coefficient
#' values \eqn{\phi_{ij}^{(\ell)}} for each variable pair as a function of \eqn{\lambda}.
#'
#' @param fit A fitted object returned by \code{\link{sparseVAR}} (or \code{cv.sparseVAR})
#'   with class \code{"bigtime.VAR"}.
#' @param Y Optional numeric matrix of time series data. If \code{NULL}, uses \code{fit$Y}.
#' @param VARpen Optional character string specifying the VAR penalty. If \code{NULL},
#'   uses \code{fit$VARpen}.
#' @param log_lambda Logical; if \code{TRUE}, the x-axis displays \eqn{\log(\lambda)}. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed down to \code{\link{sparseVAR}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing individual coefficient paths faceted by variable pair.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline labs theme_minimal facet_wrap
plot_coefficient_path <- function(fit, Y = NULL, VARpen = NULL, log_lambda = TRUE, ...) {
  
  if (is.null(Y)) {
    Y <- fit$Y
    if (is.null(Y)) {
      stop("Data matrix Y was not found in 'fit' object. Please provide Y explicitly.")
    }
  }
  
  if (is.null(VARpen)) {
    VARpen <- fit$VARpen
    if (is.null(VARpen)) {
      stop("VAR penalty was not found in 'fit' object. Please provide 'VARpen' explicitly.")
    }
  }
  
  lambda_grid <- fit$lambdas
  if (is.null(lambda_grid)) {
    stop("No lambda grid found in the provided fit object.")
  }
  
  p <- fit$p
  k <- fit$k
  series_names <- if (!is.null(fit$series_names)) fit$series_names else paste0("Y", seq_len(k))
  
  # Refit across the lambda grid and extract raw coefficients for each variable pair
  coef_list <- lapply(lambda_grid, function(l) {
    fit_l <- sparseVAR(
      Y = Y,
      p = p,
      VARpen = VARpen,
      VARlseq = l,
      selection = "none",
      ...
    )
    
    phi_mat <- fit_l$Phihat
    if (length(dim(phi_mat)) == 3) {
      phi_mat <- matrix(phi_mat[, , 1], nrow = k, ncol = k * p)
    }
    
    grid_coefs <- list()
    for (lag in seq_len(p)) {
      if (is.list(phi_mat)) {
        phi_lag <- phi_mat[[lag]]
      } else if (is.matrix(phi_mat)) {
        col_indices <- ((lag - 1) * k + 1):(lag * k)
        phi_lag <- phi_mat[, col_indices, drop = FALSE]
      } else {
        stop("Unexpected structure in fit_l$Phihat.")
      }
      
      for (i in seq_len(k)) {
        for (j in seq_len(k)) {
          pair_label <- paste0(series_names[j], " -> ", series_names[i])
          
          grid_coefs[[length(grid_coefs) + 1]] <- data.frame(
            lambda = l,
            lag = factor(paste0("Lag ", lag), levels = paste0("Lag ", seq_len(p))),
            pair = factor(pair_label),
            value = phi_lag[i, j]
          )
        }
      }
    }
    
    do.call(rbind, grid_coefs)
  })
  
  plot_data <- do.call(rbind, coef_list)
  
  # Extract CV-selected lambdas if available
  opt_val <- if (!is.null(fit$lambda_opt) && length(fit$lambda_opt) == 1 && !is.na(fit$lambda_opt)) {
    fit$lambda_opt
  } else {
    NULL
  }
  se_opt_val <- if (!is.null(fit$lambda_SEopt) && length(fit$lambda_SEopt) == 1 && !is.na(fit$lambda_SEopt)) {
    fit$lambda_SEopt
  } else {
    NULL
  }
  
  if (log_lambda) {
    plot_data$x_val <- log(plot_data$lambda)
    x_label <- expression(log(lambda))
    if (!is.null(opt_val)) opt_val <- log(opt_val)
    if (!is.null(se_opt_val)) se_opt_val <- log(se_opt_val)
  } else {
    plot_data$x_val <- plot_data$lambda
    x_label <- expression(lambda)
  }
  
  p_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = value, color = lag)) +
    ggplot2::geom_hline(yintercept = 0, color = "gray60", linetype = "solid", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ pair) +
    ggplot2::labs(
      x = x_label,
      y = expression(hat(phi)[ij]^{(l)}),
      color = "Lag",
      title = "Coefficient Paths across Penalty Grid"
    ) +
    ggplot2::theme_minimal()
  
  # Add vertical reference lines for optimal (dashed) and 1-SE (dotted) lambda
  if (!is.null(opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = opt_val, linetype = "dashed", color = "black", linewidth = 1.1)
  }
  if (!is.null(se_opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = se_opt_val, linetype = "dotted", color = "black", linewidth = 1.1)
  }
  
  return(p_obj)
}