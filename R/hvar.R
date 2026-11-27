#' Sparse Estimation of the Vector AutoRegressive (VAR) Model
#' @param Y A \eqn{T} by \eqn{k} matrix of time series. If k=1, a univariate autoregressive model is estimated.
#' @param p User-specified maximum autoregressive lag order of the VAR. Typical usage is to have the program compute its own maximum lag order based on the time series length.
#' @param h Desired forecast horizon in time-series cross-validation procedure.
#' @param VARlseq User-specified grid of values for regularization parameter corresponding to sparse penalty (for "HLagFused", this corresponds to \eqn{\lambda_1} for the HLag component). Typical usage is to have the program compute
#' its own grid. Supplying a grid of values overrides this. WARNING: use with care.
#' @param VARgran User-specified vector of granularity specifications for the penalty parameter grid:  First element specifies
#' how deep the grid should be constructed. Second element specifies how many values the grid should contain.
#' @param cvcut Proportion of observations used for model estimation in the time series cross-validation procedure. The remainder is used for forecast evaluation. Redundant if selection is not "cv".
### ORIGINAL START
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
### ORIGINAL END
### IMB START
#' @param eps a small positive numeric value giving the tolerance for convergence in the proximal gradient algorithm.
#' If NULL (default), uses 1e-3 for "HLag" and "L1" (bigtime default) and 1e-6 for "LagFused" and "HLagFused".
### IMB END
### ORIGINAL START
#' @param VARpen "HLag" (hierarchical sparse penalty) or "L1" (standard lasso penalty) penalization.
### ORIGINAL END
### IMB START
#' @param VARpen "HLag" (hierarchical sparse penalty), "L1" (standard lasso penalty), "LagFused" (row-wise lag fusion penalty across consecutive lag matrices), or "HLagFused" (combined hierarchical lag sparsity and row-wise lag fusion penalty with two tuning parameters).
#' @param VARlseq2 Optional user-specified grid of values for the second regularization parameter \eqn{\lambda_2} (LagFused component) when \code{VARpen = "HLagFused"}. Ignored for other penalties.
#' @param VARgran_l2 Optional granularity vector \code{c(depth, n_grid)} for the \eqn{\lambda_2} grid when \code{VARpen = "HLagFused"}. If NULL (default), inherits \code{VARgran}.
### IMB END
#' @param selection One of "none" (default), "cv" (Time Series Cross-Validation), "bic", "aic", "hq". Used to select the optimal penalization.
#' @param check_std Check whether data is standardised. Default is TRUE and is not recommended to be changed
#' @param verbose Logical to print value of information criteria for each lambda together with selection. Default is FALSE
#' @export
#' @return A list with the following components
#' \item{Y}{\eqn{T} by \eqn{k} matrix of time series.}
#' \item{k}{Number of time series.}
#' \item{p}{Maximum autoregressive lag order of the VAR.}
#' \item{Phihat}{Matrix of estimated autoregressive coefficients of the VAR.}
#' \item{phi0hat}{vector of VAR intercepts.}
#' \item{series_names}{names of time series}
#' \item{lambdas}{sparsity parameter grid (for "HLagFused", a two-column matrix with columns \code{lambda1} and \code{lambda2})}
#' \item{MSFEcv}{MSFE cross-validation scores for each value of the sparsity parameter in the considered grid}
#' \item{MSFEcv_all}{MSFE cross-validation full output}
#' \item{lambda_opt}{Optimal value of the sparsity parameter as selected by the time-series cross-validation procedure}
#' \item{lambda_SEopt}{Optimal value of the sparsity parameter as selected by the time-series cross-validation
#'       procedure and after applying the one-standard-error rule.  This is the value used.}
#' \item{h}{Forecast horizon h}
#' @references Nicholson William B., Wilms Ines, Bien Jacob and Matteson David S. (2020), “High-dimensional forecasting via interpretable vector autoregression”, Journal of Machine Learning Research, 21(166), 1-52.
#' @seealso \link{lagmatrix} and \link{directforecast}
#' @examples
#' data(var.example)
#' VARfit <- sparseVAR(Y = scale(Y.var)) # sparse VAR
#' ARfit <- sparseVAR(Y=scale(Y.var[,2])) # sparse AR
### ORIGINAL START
# sparseVAR <- function(Y, p=NULL, VARpen="HLag", VARlseq=NULL, VARgran=NULL,
#                       selection = c("none", "cv", "bic", "aic", "hq"),
#                       cvcut=0.9, h=1,  eps=1e-3, check_std = TRUE, verbose = FALSE){
### ORIGINAL END

### IMB START
sparseVAR <- function(Y, p=NULL, VARpen="HLag", VARlseq=NULL, VARgran=NULL,
                      selection = c("none", "cv", "bic", "aic", "hq"),
                      cvcut=0.9, h=1, eps=NULL, check_std = TRUE, verbose = FALSE,
                      VARlseq2=NULL, VARgran_l2=NULL){
  ### IMB END
  # Check Inputs
  if(!is.matrix(Y)){
    
    if(is.vector(Y) & length(Y)>1){
      Y <- matrix(Y, ncol=1)
    }else{
      stop("Y needs to be a matrix of dimension T by k")
    }
    
  }
  
  if (check_std) .check_if_standardised(Y)
  selection = match.arg(selection) 
  
  if(nrow(Y)<10){
    stop("The time series length is too small.")
  }
  
  if(!is.null(colnames(Y))){ # time series names
    series_names <- colnames(Y)
  } else {
    series_names <- NULL
  }
  
  if(h<=0){
    stop("The forecast horizon h needs to be a strictly positive integer")
  }
  
  if(!is.null(p)){
    if(p<=0){
      stop("The maximum autoregressive order of VAR needs to be a strictly positive integer")
    }
  }
  
  if((!is.vector(VARlseq) & !is.null(VARlseq)) ){ # | length(VARlseq)==1
    stop("The regularization parameter VARlseq needs to be a vector of length=>1 or NULL otherwise")
  }
  
  ### IMB START
  if(!is.null(VARlseq) && any(VARlseq < 0)){
    stop("All values in VARlseq need to be non-negative")
  }
  
  if(!is.vector(VARlseq2) & !is.null(VARlseq2)){
    stop("The regularization parameter VARlseq2 needs to be a vector of length=>1 or NULL otherwise")
  }
  
  if(!is.null(VARlseq2) && any(VARlseq2 < 0)){
    stop("All values in VARlseq2 need to be non-negative")
  }
  ### IMB END
  
  if(any((VARgran<=0)==T)){
    stop("The granularity parameters needs to be a strictly positive integers")
  }
  
  ### IMB START
  if(!is.null(VARgran_l2) && any(VARgran_l2 <= 0)){
    stop("The granularity parameters VARgran_l2 need to be strictly positive numbers")
  }
  ### IMB END
  
  if(!((cvcut<1) & (cvcut>0))){
    stop("cvcut needs to be a number between 0 and 1")
  }
  
  ### ORIGINAL START
  # if(!is.element(VARpen, c("HLag", "L1"))){
  #   stop("The type of penalization VARpen needs to be either HLag (hierarchical sparse penalization) or L1 (standard lasso penalization)")
  # }
  ### ORIGINAL END
  
  ### IMB START
  if(!is.element(VARpen, c("HLag", "L1", "LagFused", "HLagFused"))){
    stop("The type of penalization VARpen needs to be either HLag (hierarchical sparse penalization), L1 (standard lasso penalization), LagFused (time-lag fused penalization), or HLagFused (combined HLag and LagFused penalization)")
  }
  
  if(VARpen == "HLagFused" && !(selection %in% c("none", "cv"))){
    stop("For HLagFused penalization, only selection = 'none' and selection = 'cv' are currently supported.")
  }
  
  if(VARpen != "HLagFused" && (!is.null(VARlseq2) || !is.null(VARgran_l2))){
    warning("VARlseq2 and VARgran_l2 are only used when VARpen = 'HLagFused' and will be ignored.")
  }
  ### IMB END
  
  ### ORIGINAL START
  # if(!is.null(p)){
  #   if(VARpen=="HLag" & p<=1){
  #     stop("HLag penalization is only supported for p>1. Use L1 as VARpen instead")
  #   }
  # }
  ### ORIGINAL END
  
  ### IMB START
  if(!is.null(p)){
    if((VARpen %in% c("HLag", "LagFused", "HLagFused")) & p<=1){
      stop("HLag, LagFused, and HLagFused penalization are only supported for p>1. Use L1 as VARpen instead")
    }
  }
  ### IMB END
  
  ### IMB START
  # Set penalty-dependent default convergence tolerance
  if(is.null(eps)){
    eps <- if(VARpen %in% c("LagFused", "HLagFused")) 1e-6 else 1e-3
  }
  ### IMB END
  
  if(eps<=0){
    stop("The convergence tolerance parameter eps needs to be a small positive number")
  }
  
  ### IMB START
  if((VARpen %in% c("LagFused", "HLagFused")) & eps > 1e-4){
    warning("For LagFused and HLagFused penalization, a stricter tolerance (eps <= 1e-5) is recommended to ensure accurate convergence.")
  }
  ### IMB END
  
  
  if(is.null(p)){
    p <- floor(1.5*sqrt(nrow(Y)))
  }
  
  # Regularization grid
  if(is.null(VARgran)){
    VARgran1 <- 10^2
    VARgran2 <- 10
  }else{
    VARgran1 <- VARgran[1]
    VARgran2 <- VARgran[2]
  }
  
  ### IMB START
  if(is.null(VARgran_l2)){
    VARgran1_l2 <- VARgran1
    VARgran2_l2 <- VARgran2
  }else{
    VARgran1_l2 <- VARgran_l2[1]
    VARgran2_l2 <- VARgran_l2[2]
  }
  ### IMB END
  
  ### ORIGINAL START
  # if(!is.null(VARlseq)){
  #   VARgran1 <- max(VARlseq)/min(VARlseq)
  #   VARgran2 <- length(VARlseq)
  # }
  #
  # if(length(VARlseq)==1 & selection == "cv"){ # If user only provides one lambda value, don't do cross-validation
  #   selection <- "none"
  #   warning("No cross-validation is performed since only one sparsity parameter is provided.")
  # }
  ### ORIGINAL END
  
  ### IMB START
  if(!is.null(VARlseq)){
    VARgran1 <- if(min(VARlseq) > 0) max(VARlseq)/min(VARlseq) else Inf
    VARgran2 <- length(VARlseq)
  }
  
  if(!is.null(VARlseq2)){
    VARgran1_l2 <- if(min(VARlseq2) > 0) max(VARlseq2)/min(VARlseq2) else Inf
    VARgran2_l2 <- length(VARlseq2)
  }
  
  single_grid <- if(VARpen == "HLagFused"){
    length(VARlseq) == 1 && length(VARlseq2) == 1
  }else{
    length(VARlseq) == 1
  }
  
  if(single_grid & selection == "cv"){ # If user only provides one lambda combination, don't do cross-validation
    selection <- "none"
    warning("No cross-validation is performed since only one sparsity parameter combination is provided.")
  }
  ### IMB END
  
  # Set the grid of sparsity parameters
  VARdata <- HVARmodel(Y=Y, p=p, h=h)
  k <- VARdata$k # Number of time series
  
  if(selection == "cv"){ # Time series cross-validation to get optimal sparsity parameter
    ### ORIGINAL START
    # VARcv <- HVAR_cv(Y=Y, p=p, h=h, lambdaPhiseq=VARlseq, gran1=VARgran1, gran2=VARgran2, T1.cutoff=cvcut, eps=eps, type=VARpen)
    #
    # # Var estimation with selected regularization parameter
    # VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARcv$lambda_opt_oneSE, eps=eps, type=VARpen)
    ### ORIGINAL END
    
    ### IMB START
    if(VARpen == "HLagFused"){
      VARcv <- HVAR_cv(Y=Y, p=p, h=h, lambdaPhiseq=VARlseq, gran1=VARgran1, gran2=VARgran2,
                       T1.cutoff=cvcut, eps=eps, type=VARpen,
                       lambdaPhi2seq=VARlseq2, gran1_l2=VARgran1_l2, gran2_l2=VARgran2_l2)
      
      # Var estimation with selected regularization parameter pair (lambda1, lambda2)
      VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k,
                       lambdaPhi=unname(VARcv$lambda_opt_oneSE["lambda1"]),
                       lambdaPhi2=unname(VARcv$lambda_opt_oneSE["lambda2"]),
                       eps=eps, type=VARpen)
    }else{
      VARcv <- HVAR_cv(Y=Y, p=p, h=h, lambdaPhiseq=VARlseq, gran1=VARgran1, gran2=VARgran2, T1.cutoff=cvcut, eps=eps, type=VARpen)
      
      # Var estimation with selected regularization parameter
      VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARcv$lambda_opt_oneSE, eps=eps, type=VARpen)
    }
    ### IMB END
    
  }else{ # No time series cross-validation
    
    fullY <- VARdata$fullY # response matrix
    
    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }
    fullZ <- VARdata$fullZ # design matrix
    
    ### ORIGINAL START
    # Phis <- array(NA, c(k, k*p, VARgran2)) # Estimates AR coefficients for each value in the grid
    # phi0s <- array(NA, c(k, 1, VARgran2)) # Estimates of constants for each value in the grid
    #
    # # Get lambda grid if not specified by the user
    # if(is.null(VARlseq)){
    #   jj <- .lfunction3(p, k)
    #
    #   if(VARpen=="HLag"){
    #     VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ, "HVARELEM", p, k,
    #                             MN=F, alpha=1/(k+1), C=rep(1,p))
    #   }
    #
    #   if(VARpen=="L1"){
    #     VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
    #   }
    #
    # }
    #
    # for(il in 1:length(VARlseq)){ # For now in R
    #   VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARlseq[il], eps=eps, type=VARpen)
    #   Phis[,,il] <- VARmodel$Phi
    #   phi0s[,,il] <- VARmodel$phi
    # }
    ### ORIGINAL END
    
    ### IMB START
    if(VARpen == "HLagFused"){
      jj <- .lfunction3(p, k)
      
      if(is.null(VARlseq)){
        VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ, "HVARELEM", p, k,
                                MN=F, alpha=1/(k+1), C=rep(1,p))
      }
      
      if(is.null(VARlseq2)){
        VARlseq2 <- .LambdaGridE(VARgran1_l2, VARgran2_l2, jj, fullY, fullZ, "LagFused", p, k,
                                 MN=F, alpha=1/(k+1), C=rep(1,p))
      }
      
      VARlgrid <- as.matrix(expand.grid(lambda1 = VARlseq, lambda2 = VARlseq2))
      rownames(VARlgrid) <- NULL
      ngrid <- nrow(VARlgrid)
      
      Phis <- array(NA, c(k, k*p, ngrid))
      phi0s <- array(NA, c(k, 1, ngrid))
      
      for(il in 1:ngrid){
        VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k,
                         lambdaPhi=VARlgrid[il, 1], lambdaPhi2=VARlgrid[il, 2],
                         eps=eps, type=VARpen)
        Phis[,,il] <- VARmodel$Phi
        phi0s[,,il] <- VARmodel$phi
      }
      VARlout <- VARlgrid
    }else{
      # Get lambda grid if not specified by the user
      if(is.null(VARlseq)){
        jj <- .lfunction3(p, k)
        
        if(VARpen=="HLag"){
          VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ, "HVARELEM", p, k,
                                  MN=F, alpha=1/(k+1), C=rep(1,p))
        }
        
        if(VARpen == "LagFused"){
          VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ, "LagFused", p, k,
                                  MN=F, alpha=1/(k+1), C=rep(1,p))
        }
        
        if(VARpen=="L1"){
          VARlseq <- .LambdaGridE(VARgran1, VARgran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
        }
      }
      
      Phis <- array(NA, c(k, k*p, length(VARlseq))) # Estimates AR coefficients for each value in the grid
      phi0s <- array(NA, c(k, 1, length(VARlseq))) # Estimates of constants for each value in the grid
      
      for(il in 1:length(VARlseq)){ # For now in R
        VARmodel <- HVAR(fullY=VARdata$fullY, fullZ=VARdata$fullZ, p=VARdata$p, k=VARdata$k, lambdaPhi=VARlseq[il], eps=eps, type=VARpen)
        Phis[,,il] <- VARmodel$Phi
        phi0s[,,il] <- VARmodel$phi
      }
      VARlout <- VARlseq
    }
    ### IMB END
    
  }
  
  if(selection == "cv"){
    Phihat <- if (ncol(Y) == 1) matrix(VARmodel$Phi, nrow = 1) else VARmodel$Phi
    out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=Phihat, "phi0hat"=VARmodel$phi,
                "series_names"=series_names, "lambdas"=VARcv$lambda,
                "MSFEcv"=VARcv$MSFE_avg, "MSFE_all"=VARcv$MSFE_all,
                "lambda_SEopt"=VARcv$lambda_opt_oneSE,"lambda_opt"=VARcv$lambda_opt, "h"=h,
                "selection" = selection)
  }else{
    ### ORIGINAL START
    # out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=Phis, "phi0hat"=phi0s,
    #             "series_names"=series_names, "lambdas"=VARlseq,
    #             "MSFEcv"=NA, "MSFE_all"=NA,
    #             "lambda_SEopt"=NA,"lambda_opt"=NA, "h"=h,
    #             "selection" = selection)
    ### ORIGINAL END
    ### IMB START
    out <- list("k"=k, "Y"=Y, "p"=p, "Phihat"=Phis, "phi0hat"=phi0s,
                "series_names"=series_names, "lambdas"=VARlout,
                "MSFEcv"=NA, "MSFE_all"=NA,
                "lambda_SEopt"=NA,"lambda_opt"=NA, "h"=h,
                "selection" = selection)
    ### IMB END
  }
  class(out) <- "bigtime.VAR"
  
  # If user wants to use IC for selection
  if (selection %in% c("bic", "aic", "hq")) {
    out <- ic_selection(out, ic = selection, verbose = verbose)
    Phihat <- if (ncol(Y) == 1) matrix(out$Phihat, nrow = 1) else out$Phihat
    out$Phihat <- Phihat
  }
  out$VARpen <- VARpen
  out
}

HVARmodel<-function(Y, p, h=1){
  # Preliminaries
  k <- ncol(Y) # Number of Endogenous variables

  # Lagged predictor matrices
  DATAY <- embed(Y, dimension=p+h) # collect data to compute h-step ahead forecasts when selecting lambda
  fullY <- DATAY[, 1:k] # in its columns Y1, Y2, ... Yq
  fullZ <- t(DATAY[, (ncol(DATAY)-k*p+1):ncol(DATAY)])

  out<-list("fullY"=fullY, "fullZ"=fullZ, "k"=k, "p"=p)
}

### ORIGINAL START
# HVAR<-function(fullY, fullZ, p, k, lambdaPhi, eps=1e-5, type="HLag"){
### ORIGINAL END

### IMB START
HVAR<-function(fullY, fullZ, p, k, lambdaPhi, eps=NULL, type="HLag", lambdaPhi2=NULL){
  ### IMB END
  
  ### IMB START
  if(is.null(eps)){
    eps <- if(type %in% c("LagFused", "HLagFused")) 1e-6 else 1e-5
  }
  ### IMB END
  
  # HVAR estimates
  if(type=="HLag"){
    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }
    
    Phi <- HVARElemAlgcpp(array(0, dim=c(k,k*p+1,1)), fullY, fullZ, lambdaPhi, eps, p)
  }
  
  # L1 estimates
  if(type=="L1"){
    Phi <- lassoVARFistcpp(array(0,dim=c(k,k*p+1,1)), fullY, fullZ, lambdaPhi, eps, p)
  }
  
  ### IMB START
  # Time-Lag Fused estimates
  if(type=="LagFused"){
    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }
    
    Phi <- VARLagFusedcpp(array(0, dim=c(k, k*p+1, 1)), fullY, fullZ, lambdaPhi, eps, p)
  }
  
  # Combined HLag + LagFused estimates
  if(type=="HLagFused"){
    if(is.null(lambdaPhi2)){
      stop("For type = 'HLagFused', lambdaPhi2 must be provided.")
    }
    if(k==1){
      fullY <- matrix(fullY, ncol=1)
    }
    
    Phi <- VARHLagFusedcpp(array(0, dim=c(k, k*p+1, 1)), fullY, fullZ,
                           lambdaPhi, lambdaPhi2, eps, p)
  }
  ### IMB END
  
  resids <- t(t(fullY) - Phi[,,1]%*%rbind(rep(1, ncol(fullZ)), fullZ))
  Yhat <- t(Phi[,, 1]%*%rbind(rep(1, ncol(fullZ)), fullZ))
  
  # Output
  out <- list("Phi"=Phi[,-1,1], "phi0"=Phi[,1,1], "resids"=resids, "Yhat"=Yhat, "Y"=fullY)
  
}

### ORIGINAL START
# HVAR_cv<-function(Y, p, h=1, lambdaPhiseq=NULL, gran1 = 10^2, gran2=10, T1.cutoff=0.9, eps=1e-5, type="HLag"){
### ORIGINAL END

### IMB START
HVAR_cv<-function(Y, p, h=1, lambdaPhiseq=NULL, gran1 = 10^2, gran2=10, T1.cutoff=0.9, eps=NULL, type="HLag",
                  lambdaPhi2seq=NULL, gran1_l2=NULL, gran2_l2=NULL){
  ### IMB END
  
  ### IMB START
  if(is.null(eps)){
    eps <- if(type %in% c("LagFused", "HLagFused")) 1e-6 else 1e-5
  }
  
  if(is.null(gran1_l2)){
    gran1_l2 <- gran1
  }
  if(is.null(gran2_l2)){
    gran2_l2 <- gran2
  }
  ### IMB END
  
  # Get response and predictor matrix
  HVARmodelFIT <- HVARmodel(Y=Y, p=p, h=h)
  k <- HVARmodelFIT$k # Number of time series
  fullY <- HVARmodelFIT$fullY # response matrix
  if(k==1){
    fullY <- matrix(fullY, ncol=1)
  }
  fullZ <- HVARmodelFIT$fullZ # design matrix
  
  
  ### ORIGINAL START
  # # Get lambda grid if not specified by the user
  # if(is.null(lambdaPhiseq)){
  #   jj <- .lfunction3(p, k)
  #
  #   if(type=="HLag"){
  #     lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"HVARELEM", p, k,
  #                                      MN=F, alpha=1/(k+1), C=rep(1,p))
  #   }
  #
  #   if(type=="L1"){
  #     lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
  #   }
  #
  # }
  ### ORIGINAL END
  
  ### IMB START
  # Get lambda grid(s) if not specified by the user
  if(type == "HLagFused"){
    jj <- .lfunction3(p, k)
    
    if(is.null(lambdaPhiseq)){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ, "HVARELEM", p, k,
                                   MN=F, alpha=1/(k+1), C=rep(1,p))
    }
    
    if(is.null(lambdaPhi2seq)){
      lambdaPhi2seq <- .LambdaGridE(gran1_l2, gran2_l2, jj, fullY, fullZ, "LagFused", p, k,
                                    MN=F, alpha=1/(k+1), C=rep(1,p))
    }
    
    lambda_grid <- as.matrix(expand.grid(lambda1 = lambdaPhiseq, lambda2 = lambdaPhi2seq))
    rownames(lambda_grid) <- NULL
  }else if(is.null(lambdaPhiseq)){
    jj <- .lfunction3(p, k)
    
    if(type=="HLag"){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"HVARELEM", p, k,
                                   MN=F, alpha=1/(k+1), C=rep(1,p))
    }
    
    if(type == "LagFused"){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ, "LagFused", p, k,
                                   MN=F, alpha=1/(k+1), C=rep(1,p))
    }
    
    if(type=="L1"){
      lambdaPhiseq <- .LambdaGridE(gran1, gran2, jj, fullY, fullZ,"Basic",p,k,MN=F,alpha=1/(k+1),C=rep(1,p))
    }
  }
  ### IMB END
  
  # Time Series cross-validation loop
  # Choose sparsity parameter that minimizes h-step ahead mean squared prediction error.
  n <- nrow(fullY)
  T1 <- floor(T1.cutoff*n)
  tseq <- T1:(n-1)
  
  ### ORIGINAL START
  # estim <- (type=="HLag")*2 + (type=="L1")*1
  ### ORIGINAL END
  
  ### IMB START
  # Note: estim == 4 is reserved for LagFusedZero
  estim <- (type=="HLag")*2 + (type=="L1")*1 + (type=="LagFused")*3 + (type=="HLagFused")*5
  ### IMB END
  
  ### ORIGINAL START
  # # Time-series cross-validation for tuning parameter selection
  # my_cv <- HVAR_cvaux_loop_cpp(Y = fullY, Z = fullZ, tseq = tseq, gamm = lambdaPhiseq,  eps = eps, p = p, estim = estim)
  # MSFEmatrix <- my_cv$MSFEcv
  # rownames(MSFEmatrix) <- paste0("t=", tseq)
  # colnames(MSFEmatrix) <- paste0("lambda=", lambdaPhiseq)
  # sparsitymatrix <- my_cv$sparsitycv
  # colnames(sparsitymatrix) <- colnames(MSFEmatrix)
  # rownames(sparsitymatrix) <- rownames(MSFEmatrix)
  # MSFE_avg <- apply(MSFEmatrix, 2, mean)
  # lambda_opt <- lambdaPhiseq[which.min(MSFE_avg)]
  #
  # if(length(lambda_opt)==0){
  #   lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
  # }else{
  #   if(is.na(lambda_opt)){
  #     lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
  #   }
  # }
  #
  # # One-standard error rule to determine optimal lambda values
  # lambda_opt_oneSE <- NA
  # lambda_optinit <- lambda_opt
  # MSFE_sd <- apply(MSFEmatrix, 2, sd)/sqrt(nrow(MSFEmatrix))
  # MSFE_flagOK <- MSFE_avg < min(MSFE_avg, na.rm=T) + MSFE_sd[which.min(MSFE_avg)]
  # MSFEnew <- MSFE_avg
  # MSFEnew[!MSFE_flagOK] <- NA
  # sparsitynew <- apply(sparsitymatrix,2,mean)
  # sparsitynew[!MSFE_flagOK] <- NA
  #
  # lambda_opt_oneSE <- lambdaPhiseq[which.min(sparsitynew)]
  #
  # if(length(lambda_opt_oneSE)==0){
  #   lambda_opt_oneSE <- lambda_optinit
  # }else{
  #   if(is.na(lambda_opt_oneSE)){
  #     lambda_opt_oneSE <- lambda_optinit
  #   }
  # }
  #
  # gridflag_oneSE <- lambda_opt_oneSE==min(lambdaPhiseq)
  # if(lambda_opt_oneSE==min(lambdaPhiseq)){
  #   warning("Lower bound of lambda grid is selected: higher value of first granularity parameter is recommended as an input")
  # }
  #
  #
  # # Output
  # out <- list("lambda" = lambdaPhiseq,
  #             "lambda_opt_oneSE" = lambda_opt_oneSE,
  #             "MSFE_all" = MSFEmatrix,
  #             "MSFE_avg" = MSFE_avg,
  #             "flag_oneSE" = gridflag_oneSE,
  #             "lambda_opt" = lambda_opt)
  ### ORIGINAL END
  
  ### IMB START
  if(type == "HLagFused"){
    # Time-series cross-validation for 2D tuning parameter grid (lambda1, lambda2)
    my_cv <- HVAR_cvaux_loop_cpp(Y = fullY, Z = fullZ, tseq = tseq,
                                 gamm = lambda_grid[, 1], eps = eps, p = p, estim = estim,
                                 gamm2 = lambda_grid[, 2])
    MSFEmatrix <- my_cv$MSFEcv
    rownames(MSFEmatrix) <- paste0("t=", tseq)
    colnames(MSFEmatrix) <- paste0("lambda1=", lambda_grid[, 1], ",lambda2=", lambda_grid[, 2])
    sparsitymatrix <- my_cv$sparsitycv
    colnames(sparsitymatrix) <- colnames(MSFEmatrix)
    rownames(sparsitymatrix) <- rownames(MSFEmatrix)
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
    
    opt_idx <- which.min(MSFE_avg)
    if(length(opt_idx) == 0 || is.na(opt_idx)){
      opt_idx <- round(median(1:nrow(lambda_grid)))
    }
    lambda_opt <- c(lambda1 = unname(lambda_grid[opt_idx, 1]),
                    lambda2 = unname(lambda_grid[opt_idx, 2]))
    
    # One-standard error rule: sparsest model (smallest effective df) within 1 SE of minimum MSFE;
    # ties broken by lowest MSFE
    MSFE_sd <- apply(MSFEmatrix, 2, sd)/sqrt(nrow(MSFEmatrix))
    MSFE_flagOK <- MSFE_avg <= min(MSFE_avg, na.rm=T) + MSFE_sd[opt_idx]
    sparsitynew <- apply(sparsitymatrix, 2, mean)
    sparsitynew[!MSFE_flagOK] <- NA
    
    ok_idx <- which(MSFE_flagOK)
    if(length(ok_idx) > 0 && !all(is.na(sparsitynew[ok_idx]))){
      min_sp <- min(sparsitynew[ok_idx], na.rm = TRUE)
      cand_idx <- ok_idx[which(abs(sparsitynew[ok_idx] - min_sp) <= 1e-8)]
      oneSE_idx <- cand_idx[which.min(MSFE_avg[cand_idx])]
    }else{
      oneSE_idx <- opt_idx
    }
    
    lambda_opt_oneSE <- c(lambda1 = unname(lambda_grid[oneSE_idx, 1]),
                          lambda2 = unname(lambda_grid[oneSE_idx, 2]))
    
    gridflag_oneSE <- c(
      lambda1 = unname(lambda_opt_oneSE["lambda1"] == min(lambdaPhiseq)),
      lambda2 = unname(lambda_opt_oneSE["lambda2"] == min(lambdaPhi2seq))
    )
    if(gridflag_oneSE["lambda1"] && length(lambdaPhiseq) > 1 && min(lambdaPhiseq) > 0){
      warning("Lower bound of lambda1 (HLag) grid is selected: higher value of first granularity parameter is recommended as an input")
    }
    if(gridflag_oneSE["lambda2"] && length(lambdaPhi2seq) > 1 && min(lambdaPhi2seq) > 0){
      warning("Lower bound of lambda2 (LagFused) grid is selected: higher value of first granularity parameter for VARgran_l2 is recommended as an input")
    }
    
    out <- list("lambda" = lambda_grid,
                "lambda_opt_oneSE" = lambda_opt_oneSE,
                "MSFE_all" = MSFEmatrix,
                "MSFE_avg" = MSFE_avg,
                "flag_oneSE" = gridflag_oneSE,
                "lambda_opt" = lambda_opt)
  }else{
    # Time-series cross-validation for 1D tuning parameter selection
    my_cv <- HVAR_cvaux_loop_cpp(Y = fullY, Z = fullZ, tseq = tseq, gamm = lambdaPhiseq, eps = eps, p = p, estim = estim)
    MSFEmatrix <- my_cv$MSFEcv
    rownames(MSFEmatrix) <- paste0("t=", tseq)
    colnames(MSFEmatrix) <- paste0("lambda=", lambdaPhiseq)
    sparsitymatrix <- my_cv$sparsitycv
    colnames(sparsitymatrix) <- colnames(MSFEmatrix)
    rownames(sparsitymatrix) <- rownames(MSFEmatrix)
    MSFE_avg <- apply(MSFEmatrix, 2, mean)
    lambda_opt <- lambdaPhiseq[which.min(MSFE_avg)]
    
    if(length(lambda_opt)==0){
      lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
    }else{
      if(is.na(lambda_opt)){
        lambda_opt <- lambdaPhiseq[round(median(1:length(lambdaPhiseq)))]
      }
    }
    
    # One-standard error rule to determine optimal lambda values
    lambda_opt_oneSE <- NA
    lambda_optinit <- lambda_opt
    MSFE_sd <- apply(MSFEmatrix, 2, sd)/sqrt(nrow(MSFEmatrix))
    MSFE_flagOK <- MSFE_avg < min(MSFE_avg, na.rm=T) + MSFE_sd[which.min(MSFE_avg)]
    MSFEnew <- MSFE_avg
    MSFEnew[!MSFE_flagOK] <- NA
    sparsitynew <- apply(sparsitymatrix,2,mean)
    sparsitynew[!MSFE_flagOK] <- NA
    
    if (type == "LagFused") {
      # Paper 1-SE rule: select largest lambda within one standard error of minimum MSFE
      ok_idx <- which(MSFE_flagOK)
      if (length(ok_idx) > 0) {
        lambda_opt_oneSE <- max(lambdaPhiseq[ok_idx])
      } else {
        lambda_opt_oneSE <- lambda_optinit
      }
    } else {
      lambda_opt_oneSE <- lambdaPhiseq[which.min(sparsitynew)]
    }
    
    if(length(lambda_opt_oneSE)==0){
      lambda_opt_oneSE <- lambda_optinit
    }else{
      if(is.na(lambda_opt_oneSE)){
        lambda_opt_oneSE <- lambda_optinit
      }
    }
    
    gridflag_oneSE <- lambda_opt_oneSE==min(lambdaPhiseq)
    if(lambda_opt_oneSE==min(lambdaPhiseq)){
      warning("Lower bound of lambda grid is selected: higher value of first granularity parameter is recommended as an input")
    }
    
    out <- list("lambda" = lambdaPhiseq,
                "lambda_opt_oneSE" = lambda_opt_oneSE,
                "MSFE_all" = MSFEmatrix,
                "MSFE_avg" = MSFE_avg,
                "flag_oneSE" = gridflag_oneSE,
                "lambda_opt" = lambda_opt)
  }
  ### IMB END
  
}
