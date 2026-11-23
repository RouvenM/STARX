.lfunction3 <- function(p, k)
{

  kk <- .lfunction2(p,k)

  oo <- list()

  pp <- list()

  for(i in 1:length(kk))
  {
    j <- 0
    oo[[i]] <- kk[[i]][(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    pp[[i]] <- kk[[i]][-(seq(1,length(kk[[1]]),k+1)+(j*k^2))]
    j <- j+1

  }

  ownoth <- c(oo,pp)
  return(ownoth)
}


.lfunction2 <- function(p, k)
{

  kk <- list()

  kk[[1]] <- 1:(k^2)

  if(p>1)
  {
    for(i in 2:p)
    {
      kk[[i]] <- 1:(k^2)+tail(kk[[i-1]],1)

    }
  }
  return(kk)
}

### ORIGINAL START
# .LambdaGridE <- function (gran1, gran2, jj = jj, Y, Z, group, p, k, MN, alpha, C)
# {
# 
# 
#   if (group == "Basic"|group=="Tapered") {
#     gamstart <- max(t(Y) %*% t(Z))
#   }
# 
# 
#   if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM"|group=="HVARELEMSQRTLAS"|group=="HVARELEMnucl") {
#     gmax <- c()
#     for (i in 1:k) {
#       gmax[i]  <- norm2(Z %*% Y[, i])
#     }
#     gamstart <- max(gmax)
#   }
# 
#   if(group=="Tapered"){
#     beta <- array(0,dim=c(k,k*p+1,10))
#   }else{
#     beta <- array(0,dim=c(k,k*p+1,1))
#   }
# 
#   if (group == "HVARELEM" | group =="Basic"){
#     estim <- ifelse(group == "Basic", 1, 2)
#     gamstart <- LGSearch_cpp(gamstart, Y, Z, beta, estim, k, p)
#   }
#   else{
#     gamstart <- LGSearch(gamstart,Y,Z,beta,group,k,p,jj,MN,alpha,C)
#   }
# 
#   gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), length = gran2))
#   return(gamm)
# }
### ORIGINAL END

.LambdaGridE <- function (gran1, gran2, jj = jj, Y, Z, group, p, k, MN, alpha, C)
{
  if (group == "Basic"|group=="Tapered") {
    gamstart <- max(t(Y) %*% t(Z))
  }
  
  if (group == "HVARC"|group=="HVAROO"|group=="HVARELEM"|group=="HVARELEMSQRTLAS"|group=="HVARELEMnucl") {
    gmax <- c()
    for (i in 1:k) {
      gmax[i]  <- norm2(Z %*% Y[, i])
    }
    gamstart <- max(gmax)
  }
  
  if (group == "LagFused") {
    # Ensure Y is a matrix (handles k = 1) and center Y and Z as in VARLagFusedcpp
    Y_mat <- if (is.matrix(Y)) Y else matrix(Y, ncol = k)
    Y_c <- scale(Y_mat, center = TRUE, scale = FALSE)
    Z_c <- Z - rowMeans(Z)
    
    # Build backward-cumulative lag blocks Z_tilde[[m]] = sum_{l=m}^p Z_c^(l)
    Z_tilde <- vector("list", p)
    Z_tilde[[p]] <- Z_c[((p - 1) * k + 1):(p * k), , drop = FALSE]
    for (m in (p - 1):1) {
      Z_m <- Z_c[((m - 1) * k + 1):(m * k), , drop = FALSE]
      Z_tilde[[m]] <- Z_tilde[[m + 1]] + Z_m
    }
    
    # Compute pooled OLS estimates (corresponds to fully fused model where delta^(2..p) = 0)
    Z_pooled <- Z_tilde[[1]]
    Phi_bar_t <- solve(Z_pooled %*% t(Z_pooled), Z_pooled %*% Y_c)
    Res <- Y_c - t(Z_pooled) %*% Phi_bar_t
    
    # Exact KKT lambda_max: max_{i, m >= 2} || Z_tilde^(m) * r_i ||_2
    gamstart <- 0
    for (m in 2:p) {
      G_m <- Z_tilde[[m]] %*% Res
      max_norm_m <- max(sqrt(colSums(G_m^2)))
      if (max_norm_m > gamstart) {
        gamstart <- max_norm_m
      }
    }
  }

  if(group=="Tapered"){
    beta <- array(0,dim=c(k,k*p+1,10))
  }else{
    beta <- array(0,dim=c(k,k*p+1,1))
  }
 
  if (group == "HVARELEM" | group == "Basic") {
    estim <- ifelse(group == "Basic", 1, 2)
    gamstart <- LGSearch_cpp(gamstart, Y, Z, beta, estim, k, p)
  } else if (group != "LagFused") {
    gamstart <- LGSearch(gamstart, Y, Z, beta, group, k, p, jj, MN, alpha, C)
  }
  
  gamm <- exp(seq(from = log(gamstart), to = log(gamstart/gran1), length = gran2))
  return(gamm)
}
### IMB END



LGSearch<- function(gstart, Y, Z, BOLD, group, k, p, gs, MN, alpha, C)
{

  # tk <- 1/max(Mod(eigen(Z%*%t(Z))$values))
  lambdah <- gstart
  lambdal <- 0
  activeset <- list(rep(rep(list(0), length(gs))))


  while(max(abs(lambdah-lambdal))>.00001)
  {

    lambda <- (lambdah+lambdal)/2


    if(group=="Basic"){
      param <- lassoVARFistcpp(BOLD, Y, Z, lambda, 1e-04, p)[,2:(k*p+1),]
    }



    if(group=="HVARELEM")
    {
      BOLD <- HVARElemAlgcpp(BOLD, Y, Z, lambda, 1e-4, p)
      param <- BOLD[,2:(k*p+1),]
    }


    if(MN){

      diag(param[1:k,1:k]) <- ifelse(C==0,diag(param[1:k,1:k]),0)
      diag(BOLD[,2:(k*p+1),]) <- ifelse(C==0,diag(BOLD[,2:(k*p+1),]),0)

    }

    if(max(abs(param))==0)
    {
      lambdah <- lambda

    }else{
      lambdal <- lambda
    }
  }

  lambdah

}
