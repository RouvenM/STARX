#' Plot the Cross Validation Error Curve for a Sparse VAR or VARX
#' @param fit Fitted VAR, VARMA or VARX model.
#' returned by \code{\link{sparseVAR}}, \code{\link{sparseVARMA}} or \code{\link{sparseVARX}}.
#' @param ... Not currently used
#' @export
plot_cv <- function(fit, ...) {
  if ("bigtime.VAR" %in% class(fit)) {
    if (fit$selection != "cv") stop("No cross-validation was used in model estimation. Set selection='cv' in sparseVAR")
    # This copied mostly from ggb package
    m <- fit$MSFEcv
    se <- apply(fit$MSFE_all, 2, sd) / sqrt(nrow(fit$MSFE_all))
    yrang = range(c(m - se, m + se), na.rm = TRUE)
    par(xpd = FALSE)
    graphics::plot(fit$lambdas, m, xlab = "lambda",
                   ylab = "Cross-validation Error",
                   pch = 19, col = 2, ylim = yrang, log = "x", ...)
    error_bars(fit$lambdas, m - se, m + se, width = 0.01,
               col = "darkgrey")
    # graphics::abline(v = c(fit$lambda_opt, fit$lambda_SEopt), lty = 3)
    graphics::abline(v = c(fit$lambda_opt), lty = 'dashed', lwd = 2)
    graphics::abline(v = c(fit$lambda_SEopt), lty = 'dotted', lwd = 2, col = 2)
    graphics::legend("top", xpd = TRUE,
                     inset = c(0.05, -0.25),
                     legend = c("Opt. lambda", "one-SE lambda"),
                     col = c(1, 2), lty = c(3, 3), bty = "n", horiz = TRUE, lwd = c(2, 2))
  }
  else if ("bigtime.VARX" %in% class(fit)) {
    if (fit$selection != "cv") stop("No cross-validation was used in model estimation. Set selection='cv' in sparseVAR")
    image_cv(fit)
  } else if ("bigtime.VARMA" %in% class(fit)) {
    if (fit$VARMAselection != "cv") stop("No cross-validation was used during model building.")
    fit_tmp <- fit
    fit_tmp$lambdaB <- fit$PhaseII_lambdaTheta
    fit_tmp$lambdaPhi <- fit$PhaseII_lambdaPhi
    fit_tmp$MSFEcv <- fit$PhaseII_MSFEcv
    fit_tmp$lambdaPhi_SEopt <- fit$PhaseII_lambdaPhi_SEopt
    fit_tmp$lambdaB_SEopt <- fit$PhaseII_lambdaTheta_SEopt
    image_cv(fit_tmp, xlab = "lambda_Theta")
  }

  else stop("Unsupported type.")
  invisible()
}

error_bars <- function (x, upper, lower, width = 0.02, ...) {
  # this function came from hierNet package
  xlim <- range(x)
  barw <- diff(xlim) * width
  graphics::segments(x, upper, x, lower, ...)
  graphics::segments(x - barw, upper, x + barw, upper, ...)
  graphics::segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

image_cv <- function(fit, xlab = "lambda_B", ylab = "lambda_Phi") {
  o <- order(unique(fit$lambdaB))
  graphics::image(x = unique(fit$lambdaB)[o],
        y = unique(fit$lambdaPhi),
        z = matrix(fit$MSFEcv, nrow = length(o))[o, ],
        xlab = xlab,
        ylab=ylab,
        col = heat.colors(40, rev = TRUE))
  points(fit$lambdaB_SEopt, fit$lambdaPhi_SEopt, pch = 19)
}

### IMB START
#' Plot Cross-Validation Error Curve via ggplot2
#'
#' Plots the cross-validation error (MSFE) across the lambda grid with error bars (standard error)
#' in light grey tones, with high-contrast vertical lines for optimal and 1-SE optimal lambdas.
#'
#' @param fit A fitted object returned by \code{\link{sparseVAR}} (or \code{cv.sparseVAR}) 
#'   with class \code{"bigtime.VAR"}.
#' @param log_lambda Logical; if \code{TRUE}, the x-axis displays \eqn{\log(\lambda)}. Defaults to \code{TRUE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing the CV error curve.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_pointrange geom_vline labs theme_minimal
plot_cv_gg2 <- function(fit, log_lambda = TRUE, ...) {
  
  # 1. Input checks
  if (!inherits(fit, "bigtime.VAR")) {
    stop("Currently only objects of class 'bigtime.VAR' are supported by plot_cv_gg2.")
  }
  
  if (is.null(fit$selection) || fit$selection != "cv") {
    stop("No cross-validation was used in model estimation. Set selection = 'cv' in sparseVAR.")
  }
  
  lambda_grid <- fit$lambdas
  if (is.null(lambda_grid)) {
    stop("No lambda grid found in the provided fit object.")
  }
  
  # 2. Compute mean CV error (MSFE) and Standard Error
  m <- fit$MSFEcv
  se <- apply(fit$MSFE_all, 2, stats::sd) / sqrt(nrow(fit$MSFE_all))
  
  plot_data <- data.frame(
    lambda = lambda_grid,
    msfe = m,
    lower = m - se,
    upper = m + se
  )
  
  # 3. Handle x-axis scaling
  if (log_lambda) {
    plot_data$x_val <- log(plot_data$lambda)
    x_label <- expression(log(lambda))
    opt_val <- if (!is.null(fit$lambda_opt)) log(fit$lambda_opt) else NULL
    se_opt_val <- if (!is.null(fit$lambda_SEopt) && !is.na(fit$lambda_SEopt)) log(fit$lambda_SEopt) else NULL
  } else {
    plot_data$x_val <- plot_data$lambda
    x_label <- expression(lambda)
    opt_val <- fit$lambda_opt
    se_opt_val <- fit$lambda_SEopt
  }
  
  # 4. Construct plot (Light grey curve with dark red/navy lambda lines)
  p_obj <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = msfe)) +
    ggplot2::geom_line(color = "grey65", linewidth = 0.7) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = lower, ymax = upper),
      color = "grey55",
      fill = "grey55",
      size = 0.45,
      shape = 21
    ) +
    ggplot2::labs(
      x = x_label,
      y = "Cross-Validation Error (MSFE)",
      title = "Cross-Validation Error across Penalty Grid"
    ) +
    ggplot2::theme_minimal()
  
  # Optimal lambda (dunkelrot, gestrichelt)
  if (!is.null(opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = opt_val, linetype = "dashed", color = "black", linewidth = 1.1)
  }
  
  # 1-SE optimal lambda (marineblau, gepunktet)
  if (!is.null(se_opt_val)) {
    p_obj <- p_obj + ggplot2::geom_vline(xintercept = se_opt_val, linetype = "dotted", color = "black", linewidth = 1.1)
  }
  
  return(p_obj)
}
### IMB END
