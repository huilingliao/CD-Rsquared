library(Rcpp)
library(glmnet)
library(pbmcapply)
library(pbapply)
sourceCpp("matrix_multi.cpp")

## allele qc
allele.qc = function(a1,a2,ref1,ref2) {
  ref = ref1
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip1 = flip
  
  ref = ref2
  flip = ref
  flip[ref == "A"] = "T"
  flip[ref == "T"] = "A"
  flip[ref == "G"] = "C"
  flip[ref == "C"] = "G"
  flip2 = flip;
  
  snp = list()
  snp[["keep"]] = (!((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))) *
    ((a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1) | ((a1 == ref1 & a2 == ref2)) | (a1 == flip1 & a2 == flip2))
  snp[["keep"]] = as.logical(snp[["keep"]])
  snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  
  return(snp)
}


## Data simulation
#'@param n1 sample size in stage 1
#'@param r12 ratio of n2/n1, need to be >=1
#'@param p.XYU number of IVs for exposure X, outcome Y and unobserved confounder U
#'@param alpha.YtoX effect of Y on X
#'@param alpha.ZtoX effect of Z on X
#'@param alpha.UtoX effect of U on X
#'@param beta.XtoY effect of X on Y
#'@param beta.ZtoY effect of Z on Y
#'@param beta.UtoY effect of U on Y
#'@param eta.ZtoU effect of invalid IVs on U
#'@param error.X standard deviation for X
#'@param error.Y standard deviation for Y
#'@param error.U standard deviation for U
#'@param MAF.mag Minor Allele Frequence for generation of IVs
#'@param seed seed for reproduction
simu_Data <- function(n1, r12, p.XYU, alpha.YtoX, alpha.ZtoX, alpha.UtoX, beta.XtoY, beta.ZtoY, beta.UtoY, eta.ZtoU, error.X = 1, error.Y = 1, error.U = sqrt(2), MAF.mag = NULL, seed = 1025) {
  p  <- sum(p.XYU, na.rm = T)
  n2 <- n1 * r12
  n  <- n1 + n2
  
  set.seed(seed)
  # generate effect of Z to X
  if (!is.na(p.XYU[1])) {
    p.snp.X <- p.XYU[1]
    if (is.na(alpha.ZtoX)) {
      alpha.effect <- runif(p.snp.X, 0.2, 0.3) * (rbinom(p.snp.X, 1, 0.5) * 2 - 1)
    }else if(length(alpha.ZtoX) == p.snp.X) {
      alpha.effect <- alpha.ZtoX
    }else if(length(alpha.ZtoX) == 2) {
      alpha.effect <- runif(p.snp.X, alpha.ZtoX[1], alpha.ZtoX[2]) * (rbinom(p.snp.X, 1, 0.5) * 2 - 1)
    }
  }else{
    p.snp.X <- 0
    alpha.effect <- NULL
  }
  
  # generate effect of Z to Y
  if (!is.na(p.XYU[2])) {
    p.snp.Y <- p.XYU[2]
    if (is.na(beta.ZtoY)) {
      beta.effect <- runif(p.snp.Y, 0.2, 0.3) * (rbinom(p.snp.Y, 1, 0.5) * 2 - 1)
    }else if(length(beta.ZtoY) == p.snp.Y) {
      beta.effect <- alpha.ZtoX
    }else if(length(beta.ZtoY) == 2) {
      beta.effect <- runif(p.snp.Y, beta.ZtoY[1], beta.ZtoY[2]) * (rbinom(p.snp.Y, 1, 0.5) * 2 - 1)
    }
  }else{
    p.snp.Y <- 0
    beta.effect <- NULL
  }
  
  # generate effect of invalid Zs
  if (!is.na(p.XYU[3])) {
    p.snp.U <- p.XYU[3]
    if (is.na(alpha.UtoX)){
      gamma.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U,1,0.5) * 2 - 1)
    }else if(length(alpha.UtoX) == 2){
      gamma.effect <- runif(p.snp.U,alpha.UtoX[1],alpha.UtoX[2])*(rbinom(p.snp.U,1,0.5)*2-1)
    }else if(length(alpha.UtoX) == 1){
      gamma.effect <- alpha.UtoX*(rbinom(p.snp.U,1,0.5)*2-1)
    }else{
      gamma.effect <- alpha.UtoX
    }
    
    if (is.na(beta.UtoY)){
      eta.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else if(length(beta.UtoY) == 2){
      eta.effect <- runif(p.snp.U, beta.UtoY[1], beta.UtoY[2]) * (rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else if(length(beta.UtoY) == 1){
      eta.effect <- beta.UtoY*(rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else{
      eta.effect <- beta.UtoY
    }
    
    if (all(is.na(eta.ZtoU))) {
      xi.effect <- runif(p.snp.U, -0.2, 0.2)
    }else{
      xi.effect <- runif(p.snp.U, eta.ZtoU[1], eta.ZtoU[2])
    }
    p.all <- p.snp.X + p.snp.Y + p.snp.U
    MAF <- rep(ifelse(is.null(MAF.mag), 0.3, MAF.mag), p.all)
    Z <- matrix(rbinom(n * p.all, 2, MAF), ncol = p.all, byrow = T)
    Z_ref <-  matrix(rbinom(10000 * p.all, 2, MAF), ncol = p.all, byrow = T)
    U <- Z %*% matrix(c(rep(0,p.snp.X),rep(0,p.snp.Y),xi.effect), ncol = 1) + rnorm(n, 0, error.U)
  }else{
    p.snp.U <- 0
    U <- rep(0, n)
    p.all <- p.snp.X + p.snp.Y
    gamma.effect <- NULL 
    eta.effect <- NULL
    xi.effect <- NULL
    
    MAF <- rep(ifelse(is.null(MAF.mag), 0.3, MAF.mag), p.all)
    Z <- matrix(rbinom(n * p.all, 2, MAF), ncol = p.all, byrow = T)
    Z_ref <-  matrix(rbinom(10000 * p.all, 2, MAF), ncol = p.all, byrow = T)
  }
  
  eps.X <- rnorm(n, 0, error.X)
  eps.Y <- rnorm(n, 0, error.Y)
  
  X <- 0
  Y <- 0
  
  if (!is.null(alpha.effect)) {
    X <- Z[, 1:p.snp.X] %*% matrix(alpha.effect, ncol = 1)
    Y <- Z[, 1:p.snp.X] %*% matrix(beta.XtoY * alpha.effect, ncol = 1)
  }
  
  if (!is.null(beta.effect)) {
    X <- X + Z[, (p.snp.X + 1):(p.snp.X + p.snp.Y)] %*% matrix(alpha.YtoX * beta.effect, ncol = 1)
    Y <- Y + Z[, (p.snp.X + 1):(p.snp.X + p.snp.Y)] %*% matrix(beta.effect, ncol = 1)
  }
  
  if (!is.null(gamma.effect)) {
    if(!is.null(eta.effect)){
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(gamma.effect + alpha.YtoX * eta.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(beta.XtoY * gamma.effect + eta.effect, ncol = 1)
    }else{
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(gamma.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(beta.XtoY * gamma.effect, ncol = 1)
    }
  }else{
    if(!is.null(eta.effect)){
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(alpha.YtoX * eta.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(eta.effect, ncol = 1)
    }
  }
  
  X = (X + (1 + alpha.YtoX) * U + eps.X + alpha.YtoX * eps.Y)/(1 - beta.XtoY * alpha.YtoX)
  Y = (Y +  (1 + beta.XtoY) * U + beta.XtoY * eps.X + eps.Y)/(1 - beta.XtoY * alpha.YtoX)
  
  
  Z1 <- scale(Z[1:n1, ], scale = F)
  X1 <- scale(as.matrix(X[1:n1]), scale = F)
  Y1 <- scale(as.matrix(Y[1:n1]), scale = F)
  
  Z2 <- scale(Z[(n1 + 1):n, ], scale = F)
  X2 <- scale(as.matrix(X[(n1 + 1):n]), scale = F)
  Y2 <- scale(as.matrix(Y[(n1 + 1):n]), scale = F)
  
  Z_ref <- scale(Z_ref, scale = F)
  corZ_ref <- cor(Z_ref, Z_ref)
  
  ## get summary statistics
  ### for the first half
  corZ1 <- cor(Z1, Z1)
  z1Tz1 <- t(Z1) %*% Z1
  b_X <- as.numeric(1/diag(z1Tz1) * (t(Z1) %*% X1))
  rep_X1 <- X1[, rep(1, p.all)]
  se_X <- ifelse(p == 1, sqrt(colSums((X1 - Z1 * b_X)^2)/(n-2)/diag(z1Tz1)), sqrt(colSums((rep_X1 - Z1 %*% diag(b_X))^2)/(n-2)/diag(z1Tz1)))
  
  b_Y = as.numeric(1/diag(z1Tz1) * (t(Z1) %*% Y1))
  rep_Y1 = Y1[,rep(1,p.all)]
  se_Y = ifelse(p == 1,  
                sqrt(colSums((Y1 - Z1 * b_Y)^2)/(n-2)/diag(z1Tz1)), 
                sqrt(colSums((rep_Y1 - Z1 %*% diag(b_Y))^2)/(n-2)/diag(z1Tz1))
  )
  
  ### for the second half
  corZ2 <- cor(Z2, Z2)
  z2Tz2 <- t(Z2) %*% Z2
  b_X2 <- as.numeric(1/diag(z2Tz2) * (t(Z2) %*% X2))
  rep_X2 <- X2[, rep(1, p.all)]
  se_X2 <- ifelse(p == 1, sqrt(colSums((X2 - Z2 * b_X2)^2)/(n-2)/diag(z2Tz2)), sqrt(colSums((rep_X2 - Z2 %*% diag(b_X2))^2)/(n-2)/diag(z2Tz2)))
  
  b_Y2 = as.numeric(1/diag(z2Tz2) * (t(Z2) %*% Y2))
  rep_Y2 = Y2[ ,rep(1,p.all)]
  se_Y2 = ifelse(p == 1, 
                 sqrt(colSums((Y2 - Z2 * b_Y)^2)/(n-2)/diag(z2Tz2)),
                 sqrt(colSums((rep_Y2 - Z2 %*% diag(b_Y))^2)/(n-2)/diag(z2Tz2)))
  
  ## determine whether it's null or alternative 
  ## when testing X -> Y: H0: Rx2 <= Ry2 H1: Rx2 > Ry2
  ## when testing Y -> X: H0: Ry2 <= Rx2 H1: Ry2 > Rx2  
  if (alpha.YtoX == 0) {
    if (beta.XtoY == 0) {
      causal.type = "NonCausal"
    }else{
      causal.type = "XtoY"
    }
  }else{
    if (beta.XtoY == 0) {
      causal.type = "YtoX"
    }else{
      causal.type = "BiDirect"
    }
  }
  
  return(list(stage1_df = list(Y = Y1, X = X1, Z = Z1), 
              stage2_df = list(Y = Y2, X = X2, Z = Z2), 
              stat1 = list(corZ = corZ1, zTz = z1Tz1, b_X = b_X, se_X = se_X, b_Y = b_Y, se_Y = se_Y), 
              stat2 = list(corZ = corZ2, zTz = z2Tz2, b_X = b_X2, se_X = se_X2, b_Y = b_Y2, se_Y = se_Y2), 
              causal.type = causal.type, 
              corZ_ref = corZ_ref))
} 



## Simulate data from sample ADNI set
simu_ADNIdata <- function(n1, r12, cleaning_result, effect.XYU, sel_SNPs = NULL, k = 10, causal.type = "XtoY", null_Case = TRUE, seed = 1025) {
  set.seed(seed)
  n_ori <- nrow(cleaning_result)
  n <- n1 + n2
  gen_idx <- sample(n_ori, n, replace = ifelse(n > n_ori, TRUE, FALSE))
  
  beta.XY <- effect.XYU[1]
  beta.U <- effect.XYU[2]
  
  
  ## fit the original data to get coefficients
  if (is.null(sel_SNPs)) {
    sel_SNPs <- order(cor(cleaning_result[, 3], cleaning_result[, 4:33]), decreasing = T)[1:k]
  }
  
  p.all <- length(sel_SNPs)
  
  lm_formula <- paste("X ~ ", paste("SNP", sel_SNPs, sep = "",collapse =" + "))
  lm_formula <- as.formula(lm_formula)
  mdl1 <- lm(lm_formula, data = cleaning_result) # , direction = "backward", trace = -1)
  # include.idx <- sapply(strsplit(rownames(summary(lm1)$coefficients), "SNP")[-1], function(x) as.numeric(x[2]))
  
  # Predict y and u
  X_hat = predict(mdl1)
  u_hat = cleaning_result$X - X_hat
  u.sd = sd(u_hat)
  
  SNPs <- cleaning_result[gen_idx, 3 + sel_SNPs]
  SNPs <- matrix(unlist(SNPs), nrow = n)
  err.X <- rnorm(n, mean = 0, u.sd)
  err.Y <- rnorm(n, mean = 0, u.sd)
  
  ## get stage II model
  glm1 = glm(cleaning_result$Y ~ u_hat, family = binomial(link = "logit"))
  # Logistic Regression with y and u_hat
  glm2 = glm(cleaning_result$Y ~ cleaning_result$X + u_hat,
             family = binomial(link = "logit"))
  
  U <- rnorm(n, mean = 0, sd = u.sd / 2)
  X <- cbind(1, SNPs) %*% matrix(mdl1$coefficients, ncol = 1) + U * beta.U
  if (null_Case) {
    logitp <- cbind(1, U) %*% matrix(c(1, beta.U) * glm1$coefficients)
  }else{
    logitp <- cbind(1, X, U) %*% matrix(c(1, beta.XY, beta.U) * glm1$coefficients)
  }
  p <- exp(logitp) / (1 + exp(logitp))
  Y <- rbinom(n, 1, prob = p)
  
  # get datasets for both stages
  idx <- sample(n, n)
  idx1 <- idx[1:n1]
  idx2 <- idx[(n1 + 1):n]
  
  X1 <- as.matrix(X[idx1])
  X2 <- as.matrix(X[idx2])
  Y1 <- as.matrix(Y[idx1])
  Y2 <- as.matrix(Y[idx2])
  Z1 <- SNPs[idx1, ]
  Z2 <- SNPs[idx2, ]
  
  ## get summary statistics
  ### for the first half
  corZ1 <- cor(Z1, Z1)
  z1Tz1 <- t(Z1) %*% Z1
  b_X <- as.numeric(1/diag(z1Tz1) * (t(Z1) %*% X1))
  rep_X1 <- X1[, rep(1, p.all)]
  se_X <- ifelse(p.all == 1, sqrt(colSums((X1 - Z1 * b_X)^2)/(n-2)/diag(z1Tz1)), sqrt(colSums((rep_X1 - Z1 %*% diag(b_X))^2)/(n-2)/diag(z1Tz1)))
  
  b_Y = as.numeric(1/diag(z1Tz1) * (t(Z1) %*% Y1))
  rep_Y1 = Y1[,rep(1,p.all)]
  se_Y = ifelse(p.all == 1,  
                sqrt(colSums((Y1 - Z1 * b_Y)^2)/(n-2)/diag(z1Tz1)), 
                sqrt(colSums((rep_Y1 - Z1 %*% diag(b_Y))^2)/(n-2)/diag(z1Tz1))
  )
  
  ### for the second half
  corZ2 <- cor(Z2, Z2)
  z2Tz2 <- t(Z2) %*% Z2
  b_X2 <- as.numeric(1/diag(z2Tz2) * (t(Z2) %*% X2))
  rep_X2 <- X2[, rep(1, p.all)]
  se_X2 <- ifelse(p.all == 1, sqrt(colSums((X2 - Z2 * b_X2)^2)/(n-2)/diag(z2Tz2)), sqrt(colSums((rep_X2 - Z2 %*% diag(b_X2))^2)/(n-2)/diag(z2Tz2)))
  
  b_Y2 = as.numeric(1/diag(z2Tz2) * (t(Z2) %*% Y2))
  rep_Y2 = Y2[ ,rep(1,p.all)]
  se_Y2 = ifelse(p.all == 1, 
                 sqrt(colSums((Y2 - Z2 * b_Y)^2)/(n-2)/diag(z2Tz2)),
                 sqrt(colSums((rep_Y2 - Z2 %*% diag(b_Y))^2)/(n-2)/diag(z2Tz2)))
  
  ## determine whether it's null or alternative 
  ## when testing X -> Y: H0: Rx2 <= Ry2 H1: Rx2 > Ry2
  ## when testing Y -> X: H0: Ry2 <= Rx2 H1: Ry2 > Rx2  
  
  return(list(stage1_df = list(Y = Y1, X = X1, Z = Z1), 
              stage2_df = list(Y = Y2, X = X2, Z = Z2), 
              stat1 = list(corZ = corZ1, zTz = z1Tz1, b_X = b_X, se_X = se_X, b_Y = b_Y, se_Y = se_Y), 
              stat2 = list(corZ = corZ2, zTz = z2Tz2, b_X = b_X2, se_X = se_X2, b_Y = b_Y2, se_Y = se_Y2), 
              causal.type = causal.type))
}



## Get the variance of R2 in linear model
#'@param r2 Rsquared
#'@param n sample size
#'@param p number of covariates (including intercept term)
var_R2 <- function(r2, n, p) {
  # var.r2 <- 4 * r2 * (1 - r2)^2 * (n - p)^2 / ((n^2 - 1) * (n + 3))
  ## 4 / n * r2 * (1 - r2)^2 * (1 -  (2 * p + 3) / n)
  ## lbd <- n * r2 / (1 - r2)
  ## var.r2 <- 2 * ((1 - r2)^2 / n)^2 * (p + 2 * lbd)
  a <- p / 2
  b <- (n - p - 1) / 2

  # var.r2 <- ifelse(r2 > 1e-12, 4 * r2 * (1 - r2)^2 / n, a * b / ((a + b)^2 * (a + b + 1)))
  var.r2 <- ifelse(r2 > 1e-12, 4 * r2 * (1 - r2)^2 * (n - p)^2 / ((n^2 - 1) * (n + 3)), a * b / ((a + b)^2 * (a + b + 1)))
  # var.r2 <- 4 * r2 * (1 - r2)^2 / n
  return(var.r2)
}

## Testing depends on F
#'@param rx2 scalar, the R-squared for stage 1 regression
#'@param ry2 scalar, the R-squared for stage 2 regression
#'@param M the number of SNPs
#'@param n1 the sample size for stage 1
#'@param n2 the sample size for stage 2
#'@param lbd1 noncentrality parameter for stage 1 regression
#'@param lbd2 noncentrality parameter for stage 2 regression
#'@param num.sample number of samples for bootstrap
fTest <- function(rx2, ry2, M1, M2, n1, n2, num.sample = 10000) {
  statx <- rx2 * (n1 - M1 - 1) / ((1 - rx2) * M1)
  staty <- ry2 * (n2 - M2 - 1) / ((1 - ry2) * M2)
  nu1 <- max(0, (statx * M1 - M1) / n1)
  nu2 <- max(0, (staty * M2 - M2) / n2)
  # nu <- max(0, 0.5 * (lbd1 / n1 + lbd2 / n2))
  lbdx <- nu1 * n1 # nu1 * n1
  lbdy <- nu2 * n2 # nu2 * n2
  
  Tx <- rf(num.sample, M1, n1 - M1 - 1, lbdx)
  Ty <- rf(num.sample, M2, n2 - M2 - 1, lbdy)
  
  TnuX <- (Tx * M1 - M1) / n1
  TnuY <- (Ty * M2 - M2) / n2
  
  pval_xy <- mean(TnuX < TnuY)
  # CIs <- quantile(TnuX - TnuY, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  # pval_xy <- mean(statx - staty < Tx - Ty)
  # pval_yx <- mean(staty - statx < Ty - Tx)
  return(list(test.stat = statx - staty, pval = pval_xy)) # , CI = CIs))
}

## Testing depends on F approximating normal when the numbers of sample size and number of IVs go to sufficiently large
#'@param rx2 scalar, the R-squared for stage 1 regression
#'@param ry2 scalar, the R-squared for stage 2 regression
#'@param M the number of SNPs
#'@param n1 the sample size for stage 1
#'@param n2 the sample size for stage 2
#'@param lbd1 noncentrality parameter for stage 1 regression
#'@param lbd2 noncentrality parameter for stage 2 regression
#'@param num.sample number of samples for bootstrap
fTest2 <- function(rx2, ry2, M, n1, n2, lbd1, lbd2) {
  ## get the statistics
  statx <- rx2 / (1 - rx2)
  staty <- ry2 / (1 - ry2)
  
  # cat("Rx2 = ", rx2, " | Ry2 = ", ry2, "\n")
  # cat("lbd1 = ", lbd1, " | lbd2 = ", lbd2, "\n")
  varx <- 2 / (n1 - M1 - 3) / (n1 - M1 - 5) * ((M1 + lbd1)^2 / (n1 - M1 - 3) + M1 + 2 * lbd1)
  vary <- 2 / (n2 - M2 - 3) / (n2 - M2 - 5) * ((M2 + lbd2)^2 / (n2 - M2 - 3) + M2 + 2 * lbd2)
  
  # cat("var1 = ", varx, " | var2 = ", vary, "\n")
  estsd <- sqrt(varx + vary)
  
  diff_stat <- (statx - staty) / estsd
  pval_xy <- 1 - pnorm(diff_stat)
  # CIs <- c(statx - staty + qnorm(alpha/2) * estsd, statx - staty + qnorm(1 - alpha/2) * estsd)
  
  return(list(test.stat = diff_stat, pval = pval_xy)) # , CI = CIs))
}


## Calculate the ss-based and r2-based R-squared for GLM
glmR2_func <- function(y, phat) {
  SST <- sum((y - mean(y))^2)
  SSE <- sum((y - phat)^2)
  Rss2 <- 1 - SSE / SST
  Rr2 <- sum((y - mean(y)) * (phat - mean(y)))^2 / (sum((y - mean(y))^2) * sum((phat - mean(y))^2))
  return(c(Rss2, Rr2))
}

## Test based on R-squared
#'@param rx2 r-squared for exposure ~ IVs
#'@param ry2 r-squared for outcome ~ IVs
#'@param M number of IVs
#'@param n1 sample size for stage 1
#'@param n2 sample size for stage 2
#'@param test.details list of test details, list(type = "LM" or "Logistic", par = NULL or list(y, yhat))
#'@param alpha significance level of interest
R2Test <- function(rx2, ry2, M1, M2, n1, n2) {
  
  # if(type == "LM") {
    # var1 <- 4 * rx2 * (1 - rx2)^2 * (n1 - M - 1)^2 / ((n1^2 - 1) * (n1 + 3))
    # var2 <- 4 * ry2 * (1 - ry2)^2 * (n2 - M - 1)^2 / ((n2^2 - 1) * (n2 + 3))
    # var1 <- 4 / n1 * rx2 * (1 - rx2)^2 * (1 -  (2 * M + 5) / n1)
    # var2 <- 4 / n2 * ry2 * (1 - ry2)^2 * (1 -  (2 * M + 5) / n2)
    # var2 <- 4 / n2 * ry2 * (1 - ry2)^2 * (1 -  (2 * M + 5) / n2)
    # lbd1 <- n1 * rx2 / (1 - rx2)
    # lbd2 <- n2 * ry2 / (1 - ry2)
    # var1 <- 2 * ((1 - rx2)^2 / n1)^2 * (M + 2 * lbd1)
    # var2 <- 2 * ((1 - ry2)^2 / n2)^2 * (M + 2 * lbd2)
    var1 <- var_R2(rx2, n1, M1)
    var2 <- var_R2(ry2, n2, M2)
    
    # var1 <- asymR2_variance(SNPs1, rho1)
    # var2 <- asymR2_variance(SNPs2, rho2)
  # }else{
    # y <- test.details$y
    # yhat <- test.details$yhat
    
    # pp <- mean(yhat)
    # V1 <- mean(yhat * (1 - yhat))
    # V2 <- var(y)
    ## V2 <- pp * (1 - pp)
    # Z <- cbind(y, y * yhat, yhat^2)
    # Sig <- cov(Z)
    # c1 <- matrix(c((V2 - V1 + 2 * pp * V1) / V2^2, -2 / V2, 1 / V2), nrow = 1)
    # var1 <- var_R2(rx2, n1, M + 1)
    # var2 <- sqrt((c1 %*% Sig %*% t(c1))[1])
  # }
  estsd <- sqrt(var1 + var2)
  zstat <- (rx2 - ry2) / estsd
  pval_xy <- 1 - pnorm(zstat)
  
  # CIs <- c(rx2 - ry2 + qnorm(alpha/2) * estsd, rx2 - ry2 + qnorm(1 - alpha/2) * estsd)
  return(list(test.stat = rx2 - ry2, pval = pval_xy, se = estsd))
}


## Test based on R
#'@param rx correlation between X and Xhat
#'@param ry correlation between y and yhat
RTest <- function(rx, ry, n1, n2) {
  var_x <- (1 - rx^2)^2 / (n1 - 2)
  var_y <- (1 - ry^2)^2 / (n2 - 2)
  
  estsd <- sqrt(var_x + var_y)
  zstat <- (rx - ry) / estsd
  pval_xy <- 1 - pnorm(zstat)
  # CIs <- c(rx - ry + qnorm(alpha/2) * estsd, rx - ry + qnorm(1 - alpha/2) * estsd)
  return(list(test.stat = rx - ry, pval = pval_xy, se = estsd))
}

## Steiger's test
#'@param r1 correlation between exposure and IVs
#'@param r2 correlation between outcome and IVs
#'@param n1 number of sample size of stage 1
#'@param n2 number of sample size of stage 2
#'@param alpha significance level of interest
Steiger <- function(r1, r2, n1, n2) {
  p <- length(r1) 
  
  # select the SNP for test statistic construction
  r12 <- abs(r1) + abs(r2)
  idx <- which.max(r12)
  
  # construct the test statistics
  z1 <- 0.5 * log((1 + abs(r1[idx]))/(1 - abs(r1[idx])))
  z2 <- 0.5 * log((1 + abs(r2[idx]))/(1 - abs(r2[idx])))
  estsd <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z <- (z1 - z2) / estsd
  pval_xy <- 1 - pnorm(z)
  
  # CIs <- c(z1 - z2 + qnorm(alpha/2) * estsd, z1 - z2 + qnorm(1 - alpha/2) * estsd)
  return(list(test.stat = z1 - z2, pval = pval_xy, se = estsd))
}


Steiger_bonferroni <- function(r1, r2, n1, n2){
    p <- length(r1)

    z1 <- 0.5 * log((1 + abs(r1)) / (1 - abs(r1)))
    z2 <- 0.5 * log((1 + abs(r2)) / (1 - abs(r2)))

    estsd <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z <- (z1 - z2) / estsd
    cutoff_alpha <- 0.05 / p
    pval_xy <- 1 - pnorm(z)

    CIs <- cbind(z1 - z2 + qnorm(alpha/p/2) * estsd, z1 - z2 + qnorm(1 - alpha/p/2) * estsd)
    return(list(test.stat = z1 - z2, pval = pval_xy, se = estsd, CI = CIs))   
}

# CD-Ratio
#' Generate matrix V
#'
#' Generate matrix V in formula 3.6 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param SNP n by p matrix, genotype data for n individuals of p SNPs. And it is
#' standardized with column means 0 and variances 1.
#' @param rho (n+1) by (n+1) matrix, sample correlation matrix of n SNPs and X.
#'
#' @return V V matrix.
#' @export
generate_V <- function(SNP,rho)
{
  n = ncol(SNP)
  Sigma = rho[1:n,1:n]
  inv_Sigma = solve(Sigma,tol = 0)
  rho_X = as.matrix(rho[1:n,(n+1)])
  alpha =  inv_Sigma %*% rho_X
  e2 = 1 - max(0, min(1, t(rho_X) %*% inv_Sigma %*% rho_X))
  
  if(nrow(alpha)>1)
  {
    SNP_alpha = SNP %*% diag(as.numeric(alpha))
  } else{
    SNP_alpha = SNP * as.numeric(alpha)
  }
  S = rowSums(SNP_alpha)
  ###
  
  sampleszie = length(S)
  sim_X = S + rnorm(sampleszie,0,sqrt(e2))
  sim_X = scale(sim_X)*sqrt(sampleszie) / sqrt(sampleszie-1)
  M_SX = cbind(SNP,sim_X)
  
  BigM = NULL
  for(i in 1:(n+1))
  {
    BigM = cbind(BigM,M_SX[,i]*M_SX)
  }
  ###
  B = t(BigM)
  C = eigenMapMatMult(B,BigM)
  V = C / sampleszie
  ###
  #V = t(BigM)%*%BigM / sampleszie
  ###
  
  vSigma = as.matrix(c(rho))
  V = V - vSigma%*%t(vSigma)
  
  return(V)
}

#' Generate matrix M_s
#'
#' Generate matrix M_s in formula 2.9 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param n dimension.
#'
#' @return M_s M_s matrix.
#' @export
generate_M_s <- function(n)
{
  # For n, generate matrix M_s in formula 2.9
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      ind_row = (i-1)*n + j
      ind_col = (j-1)*n + i
      K[ind_row,ind_col] = K[ind_row,ind_col] + 1
    }
  }
  M_s = 1/2 * (diag(n^2) + K)
  return(M_s)
}


#' Generate matrix M_d
#'
#' Generate matrix M_d in formula 2.13 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param n dimension.
#'
#' @return M_d M_d matrix.
#' @export
generate_M_d <- function(n)
{
  # For n, generate matrix M_d in formula 2.13
  K = matrix(0,n^2,n^2)
  for(i in 1:n)
  {
    ind_rc = (i-1)*n+i
    K[ind_rc,ind_rc] = K[ind_rc,ind_rc] + 1
  }
  
  return(K)
}


#' Calculate Asymptotic Variance Matrix
#'
#' Calculate the asymptotic covariance matrix of sample correlations of n
#' SNPs with X, Use Theorem 2 in
#' "The Asymptotic Variance Matrix of the Sample Correlation Matrix".
#'
#' @param SNP n by p matrix, genotype data for n individuals of p SNPs.
#' @param rho (n+1) by (n+1) matrix, sample correlation matrix of n SNPs and X.
#'
#' @return n by n matrix, which is the asymptotic covariance matrix of
#' sample correlations of n SNPs with X.
#' @export
calculate_asymptotic_variance <- function(SNP,rho)
{
  n = ncol(rho)-1
  SNP = scale(SNP)
  SNP = SNP * sqrt(nrow(SNP)) / sqrt(nrow(SNP)-1)
  M_s = generate_M_s(n + 1)
  M_d = generate_M_d(n + 1)
  ###
  B = ( kronecker( diag(n+1) , rho) )
  C = eigenMapMatMult(M_s, B)
  C = eigenMapMatMult(C,M_d)
  M_1 = diag((n+1)^2) - C
  
  V = generate_V(SNP,rho)
  
  C = eigenMapMatMult(M_1,V)
  B = t(M_1)
  asymp_cov = eigenMapMatMult(C,B)
  
  target_ind1 = n * (n+1) + 1
  target_ind2 = (n+1)^2-1
  
  return( asymp_cov[(target_ind1:target_ind2) , (target_ind1:target_ind2)] )
}


CDratio <- function(Z, rhoxz, rhoyz, n1, n2) {
  p <- dim(Z)[2]
  
  rhox <- matrix(0, ncol = p + 1, nrow = p + 1)
  rhoy <- matrix(0, ncol = p + 1, nrow = p + 1)
  corZ <- cor(Z)
  rhox[1:p, 1:p] <- corZ
  rhox[1:p, p + 1] <- rhoxz
  rhox[p + 1, 1:p] <- rhoxz
  rhox[p + 1, p + 1] <- 1
  
  rhoy[1:p, 1:p] <- corZ
  rhoy[1:p, p + 1] <- rhoyz
  rhoy[p + 1, 1:p] <- rhoyz
  rhoy[p + 1, p + 1] <- 1
  
  V1.asy <- calculate_asymptotic_variance(Z, rhox)
  V2.asy <- calculate_asymptotic_variance(Z, rhoy)
  
  # cat("rhoyz = ", rhoyz, "\n")
  # cat("rhoxz = ", rhoxz, "\n")
  # cat("V1.asy = ", V1.asy, "\n")
  # cat("V2.asy = ", V2.asy, "\n")
  
  Kyx <- rhoyz / rhoxz
  Kxy <- rhoxz / rhoyz
  
  ## K.YX
  jc1 <- matrix(0, p, p)
  jc2 <- matrix(0, p, p)
  diag(jc1) <- 1 / rhoxz
  diag(jc2) <- -rhoyz / rhoxz^2
  jacobian <- cbind(jc1, jc2)
  combined_V <- rbind(cbind(V2.asy, matrix(0, p, p)) / n2, 
                      cbind(matrix(0, p, p), V1.asy) / n1)
  # print(jacobian)
  # print(combined_V)
  V <- jacobian %*% combined_V %*% t(jacobian)
  
  inv_V <- tryCatch(solve(V, tol = 0), error = function(e) "error")
  if (is.character(inv_V)) {
    svd_V1 <- svd(V)
    inv_V <- svd_V1$v %*% diag(1 / svd_V1$d) %*% t(svd_V1$u)
  }
  
  K.YX <- sum(inv_V %*% matrix(Kyx, ncol = 1)) / sum(inv_V)
  var.K.YX <- 1 / sum(inv_V)
  
  ## K.XY
  jc1 <- matrix(0, p, p)
  jc2 <- matrix(0, p, p)
  diag(jc1) <- 1 / rhoyz
  diag(jc2) <- - rhoxz / rhoyz^2
  jacobian <- cbind(jc1, jc2)
  combined_V <- rbind(cbind(V1.asy, matrix(0, p, p)) / n1, 
                      cbind(matrix(0, p, p), V2.asy) / n2)
  # print(jacobian)
  # print(combined_V)
  V <- jacobian %*% combined_V %*% t(jacobian)
  inv_V <- tryCatch(solve(V, tol = 0), error = function(e) "error")
  if (is.character(inv_V)) {
    svd_V1 <- svd(V)
    inv_V <- svd_V1$v %*% diag(1 / svd_V1$d) %*% t(svd_V1$u)
  }
  
  K.XY <- sum(inv_V %*% matrix(Kxy, ncol = 1)) / sum(inv_V)
  var.K.XY <- 1 / sum(inv_V)
  
  return(list(test.stat = c(K.YX, K.XY), 
              se.K.YX = sqrt(var.K.YX), 
              se.K.XY = sqrt(var.K.XY) 
  	)
  )
}


asymR2_variance <- function(SNP, rho){
  
  corZ <- cor(SNP)
  p <- length(rho)
  rhox <- matrix(NA, nrow = p + 1, ncol = p + 1)
  rhox[1:p, 1:p] <- corZ
  rhox[1:p, p + 1] <- rho
  rhox[p + 1, 1:p] <- rho
  rhox[p + 1, p + 1] <- 1
  
  V.asym <- calculate_asymptotic_variance(SNP, rhox)
  
  invCorZ <- solve(corZ, tol = 0)
  partialrho <- 2 * invCorZ %*% matrix(rho, ncol = 1)
  res <- t(partialrho) %*% V.asym %*% partialrho
  return(res)
}

genCor <- function(p, mag = 0) {
  if (mag == 0) {
    res <- diag(p)
  }else{
    res <- matrix(NA, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        res[i, j] <- exp(- abs(i - j)) * mag
      }
    }
  }
  return(res)
}

hypergeo_analytical <- function(c, z) {
  if (z == 0) {
    return(1)
  } else if (z == 1) {
    return((c - 1) / (c - 2))
  } else if (c %% 1 == 0) {
    shared <- ((z - 1) / z)
    the_sum <- 0
    for (k in seq(from = 2, to = (c - 1), by = 1)) {
      the_sum <- the_sum + shared^k / (c - k)
    }
    prefactor <- (c - 1) * z * (z - 1)^(-2)
    result <- prefactor * (the_sum - shared^c * log(1 - z))
  } else {
    cur_res <- 1 / (1 - z) * (1 + (sqrt(z) * asin(sqrt(z))) / sqrt(1 - z))
    for (i in seq(from = 1.5, to = c, by = 1)) {
      cur_c <- i - 1
      cur_res <- (cur_c - cur_c * (1 - z) * cur_res) / (z * (cur_c - 1))
    }
    result <- cur_res
  }
  return(result)
}

hypergeo_analytical_taylor <- function(c, z, tol = 1e-8, maxIter = 10000) {
  c0 <- 1
  s0 <- 1
  flag <- FALSE
  j <- 0
  
  while(!flag) {
    c1 <- c0 * (1 + j) * (1 + j) / (c + j) * z / (j + 1)
    s1 <- s0 + c1
    flag <- (abs(c1) / abs(s0) < tol | j > maxIter)
    # flag2 <- two successive terms are small compared to SN
    
    j <- j + 1
    c0 <- c1
    s0 <- s1
  }
  return(s0)
}

adjR2 <- function(X, Y, sigma_e2, R2) {
  term1 <- sum(diag(X %*% solve(t(X) %*% X, tol = 0) %*% t(X))) / (length(Y) - 1)
  res <- R2 - term1 * sigma_e2 / var(Y)
  return(res)
}


loop_func <- function(n1, r12, test.type, simuData.attrs, alpha = 0.05, MAF.mag = NULL, seed = 1025) {
  # get sample size
  n2 <- n1 * r12
  n <- n1 + n2
  
  
  ## generate data
  simuData.type = simuData.attrs$type
  if (simuData.type == "Bi") {
    p.XYU <- simuData.attrs$p.XYU
    coef.multi.X <- simuData.attrs$coef.multi.X
    coef.multi.Y <- simuData.attrs$coef.multi.Y
    multi.error <- simuData.attrs$multi.error
    mdl2.type <- simuData.attrs$mdl2.type
    p.all <- sum(p.XYU, na.rm = T)
    
    # get coefficients
    alpha.YtoX <- coef.multi.X[[1]]
    alpha.ZtoX <- coef.multi.X[[2]]
    alpha.UtoX <- coef.multi.X[[3]]
    
    beta.XtoY <- coef.multi.Y[[1]]
    beta.ZtoY <- coef.multi.Y[[2]]
    beta.UtoY <- coef.multi.Y[[3]]
    
    error.X <- multi.error[1]
    error.Y <- multi.error[2]
    error.U <- multi.error[3]
    simuData <- simu_Data(n1, r12, p.XYU, alpha.YtoX, alpha.ZtoX, alpha.UtoX, beta.XtoY, beta.ZtoY, beta.UtoY, eta.ZtoU, error.X, error.Y, error.U, MAF.mag, seed)
  }else if (simuData.type == "ADNI") {
    clearn_result <- simuData.attrs$clean_result
    effect.XYU <- simuData.attrs$effect.XYU
    sel_SNPs <- simuData.attrs$sel_SNPs
    mdl2.type <- simuData.attrs$mdl2.type
    k <- simuData.attrs$k
    causal.type <- simuData.attrs$causal.type
    simuData <- simu_ADNIdata(n1, r12, clean_result, effect.XYU, sel_SNPs, k, mdl2.type, causal.type, seed)
    p.all <- dim(simuData$stage1_df$Z)[2]
  }else if(simuData.type == "Real") {
    simuData <- simuData.attrs$Data
    p.all <- dim(simuData$stage1_df$Z)[2]
  }
  
  
  stage1_df <- simuData$stage1_df
  stage2_df <- simuData$stage2_df
  stat1 <- simuData$stat1
  stat2 <- simuData$stat2
  
  ## X to Y
  X1 <- stage1_df$X
  Y1 <- stage1_df$Y
  Z1 <- stage1_df$Z
  
  X2 <- stage2_df$X
  Y2 <- stage2_df$Y
  Z2 <- stage2_df$Z
  
  Rx2.XtoY <- summary(lm(X1 ~ Z1))$r.squared
  Ry2.ind.XtoY <- summary(lm(Y2 ~ Z2))$r.squared
  ryz <- stat2$b_Y / sqrt(stat2$b_Y^2 + (n2-1) * stat2$se_Y^2)
  CorMatInv1 <- solve(stat1$corZ, tol = 0)
  # CorMatInv1 <- tryCatch(solve(stat1$corZ), error = function(e) "error")
  # if (is.character(CorMatInv1)) {
  #   svd_Z1 <- svd(stat1$corZ)
  #   CorMatInv1 <- svd_Z1$v %*% diag(1 / svd_Z1$d) %*% t(svd_Z1$u)
  # }
  
  Ry2.sum.XtoY <- matrix(ryz, nrow = 1) %*% CorMatInv1 %*% matrix(ryz, ncol = 1)
  
  var1.XtoY <- var_R2(Rx2.XtoY, n1, p.all + 1) # 4 / n1 * Rx2.XtoY * (1 - Rx2.XtoY)^2 * (1 -  (2 * p.all + 5) / n1)
  var2.ind.XtoY <- var_R2(Ry2.ind.XtoY, n2, p.all + 1) # 4 / n2 * Ry2.ind.XtoY * (1 - Ry2.ind.XtoY)^2 * (1 -  (2 * p.all + 5) / n)
  var2.sum.XtoY <- var_R2(Ry2.sum.XtoY, n2, p.all + 1) # 4 / n2 * Ry2.sum.XtoY * (1 - Ry2.sum.XtoY)^2 * (1 -  (2 * p.all + 5) / n)
  estsd.ind.XtoY <- sqrt(var1.XtoY + var2.ind.XtoY)
  estsd.sum.XtoY <- sqrt(var1.XtoY + var2.sum.XtoY)
  
  ## Y to X
  Ry2.YtoX <- summary(lm(Y1 ~ Z1))$r.squared
  Rx2.ind.YtoX <- summary(lm(X2 ~ Z2))$r.squared
  rxz <- stat2$b_X / sqrt(stat2$b_X^2 + (n2-1) * stat2$se_X^2)
  Rx2.sum.YtoX <- matrix(rxz, nrow = 1) %*% CorMatInv1 %*% matrix(rxz, ncol = 1)
  
  var1.YtoX <- var_R2(Ry2.YtoX, n1, p.all + 1) # 4 / n1 * Rx2.XtoY * (1 - Rx2.XtoY)^2 * (1 -  (2 * p.all + 5) / n1)
  var2.ind.YtoX <- var_R2(Rx2.ind.YtoX, n2, p.all + 1) # 4 / n2 * Ry2.ind.XtoY * (1 - Ry2.ind.XtoY)^2 * (1 -  (2 * p.all + 5) / n)
  var2.sum.YtoX <- var_R2(Rx2.sum.YtoX, n2, p.all + 1) # 4 / n2 * Ry2.sum.XtoY * (1 - Ry2.sum.XtoY)^2 * (1 -  (2 * p.all + 5) / n)
  estsd.ind.YtoX <- sqrt(var1.YtoX + var2.ind.YtoX)
  estsd.sum.YtoX <- sqrt(var1.YtoX + var2.sum.YtoX)
  
  # determine H0 and H1 
  causal.type <- simuData$causal.type
  H0.XtoY.type <- c("NonCausal", "YtoX")
  H1.XtoY.type <- c("XtoY", "BiDirect")
  H0.YtoX.type <- c("NonCausal", "XtoY")
  H1.YtoX.type <- c("YtoX", "BiDirect")
  
  H0.idx.XtoY <- (causal.type %in% H0.XtoY.type)
  H0.idx.YtoX <- (causal.type %in% H0.YtoX.type)
  
  ## Hypothesis testing
  # F
  if ("F" %in% test.type) {
    # X to Y
    lbd1 <- max(0, n1 * Rx2.XtoY / (1 - Rx2.XtoY))
    lbd2.ind <- max(n2 * Ry2.ind.XtoY / (1 - Ry2.ind.XtoY))
    lbd2.sum <- max(n2 * Ry2.sum.XtoY / (1 - Ry2.sum.XtoY))
    F.res.XtoY.ind <- fTest(Rx2.XtoY, Ry2.ind.XtoY, p.all, n1, n2, alpha = alpha)
    F.res.XtoY.sum <- fTest(Rx2.XtoY, Ry2.sum.XtoY, p.all, n1, n2, alpha = alpha)
    
    if(H0.idx.XtoY) {
      type.I.err.idx.ind <- (F.res.XtoY.ind$CI[1] > 0)
      type.I.err.idx.sum <- (F.res.XtoY.sum$CI[1] > 0)
      
      F.res.XtoY <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (F.res.XtoY.ind$CI[1] > 0)
      power.idx.sum <- (F.res.XtoY.sum$CI[1] > 0)
      
      F.res.XtoY <- c(power.idx.ind, power.idx.sum)
    }
    
    # Y to X
    lbd1 <- max(0, n1 * Ry2.YtoX / (1 - Ry2.YtoX))
    lbd2.ind <- max(0, n2 * Rx2.ind.YtoX / (1 - Rx2.ind.YtoX))
    lbd2.sum <- max(0, n2 * Rx2.sum.YtoX / (1 - Rx2.sum.YtoX))
    F.res.YtoX.ind <- fTest(Ry2.YtoX, Rx2.ind.YtoX, p.all, n1, n2, alpha = alpha)
    F.res.YtoX.sum <- fTest(Ry2.YtoX, Rx2.sum.YtoX, p.all, n1, n2, alpha = alpha)
    
    if(H0.idx.YtoX) {
      type.I.err.idx.ind <- (F.res.YtoX.ind$CI[1] > 0)
      type.I.err.idx.sum <- (F.res.YtoX.sum$CI[1] > 0)
      
      F.res.YtoX <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (F.res.YtoX.ind$CI[1] > 0)
      power.idx.sum <- (F.res.YtoX.sum$CI[1] > 0)
      
      F.res.YtoX <- c(power.idx.ind, power.idx.sum)
    }
  }else{
    F.res.XtoY <- NULL
    F.res.YtoX <- NULL
    F.res.XtoY.ind <- NULL
    F.res.XtoY.sum <- NULL
    F.res.YtoX.ind <- NULL
    F.res.YtoX.sum <- NULL
  }
  
  # F2
  if ("F2" %in% test.type) {
    ## X to Y
    lbd1 <- max(0, n1 * Rx2.XtoY / (1 - Rx2.XtoY))
    lbd2.ind <- max(0, n2 * Ry2.ind.XtoY / (1 - Ry2.ind.XtoY))
    lbd2.sum <- max(0, n2 * Ry2.sum.XtoY / (1 - Ry2.sum.XtoY))
    F2.res.XtoY.ind <- fTest2(Rx2.XtoY, Ry2.ind.XtoY, p.all, n1, n2, lbd1, lbd2.ind, alpha)
    F2.res.XtoY.sum <- fTest2(Rx2.XtoY, Ry2.sum.XtoY, p.all, n1, n2, lbd1, lbd2.sum, alpha)
    
    if(H0.idx.XtoY) {
      type.I.err.idx.ind <- (F2.res.XtoY.ind$CI[1] > 0)
      type.I.err.idx.sum <- (F2.res.XtoY.sum$CI[1] > 0)
      
      F2.res.XtoY <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (F2.res.XtoY.ind$CI[1] > 0)
      power.idx.sum <- (F2.res.XtoY.sum$CI[1] > 0)
      
      F2.res.XtoY <- c(power.idx.ind, power.idx.sum)
    }
    
    ## Y to X
    lbd1 <- max(0, n1 * Ry2.YtoX / (1 - Ry2.YtoX))
    lbd2.ind <- max(0, n2 * Rx2.ind.YtoX / (1 - Rx2.ind.YtoX))
    lbd2.sum <- max(0, n2 * Rx2.sum.YtoX / (1 - Rx2.sum.YtoX))
    F2.res.YtoX.ind <- fTest2(Ry2.YtoX, Rx2.ind.YtoX, p.all, n1, n2, lbd1, lbd2.ind, alpha)
    F2.res.YtoX.sum <- fTest2(Ry2.YtoX, Rx2.sum.YtoX, p.all, n1, n2, lbd1, lbd2.sum, alpha)
    
    if(H0.idx.YtoX) {
      type.I.err.idx.ind <- (F2.res.YtoX.ind$CI[1] > 0)
      type.I.err.idx.sum <- (F2.res.YtoX.sum$CI[1] > 0)
      
      F2.res.YtoX <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (F2.res.YtoX.ind$CI[1] > 0)
      power.idx.sum <- (F2.res.YtoX.sum$CI[1] > 0)
      
      F2.res.YtoX <- c(power.idx.ind, power.idx.sum)
    }
  }else{
    F2.res.XtoY <- NULL
    F2.res.YtoX <- NULL
    F2.res.XtoY.ind <- NULL
    F2.res.XtoY.sum <- NULL
    F2.res.YtoX.ind <- NULL
    F2.res.YtoX.sum <- NULL
  }
  
  # R2
  if ("R2-LM" %in% test.type) {
    ## linear model
    test.details <- list(type = "LM")
    # X to Y
    R2.LM.res.XtoY.ind <- R2Test(Rx2.XtoY, Ry2.ind.XtoY, p.all, n1, n2, test.details, alpha)
    R2.LM.res.XtoY.sum <- R2Test(Rx2.XtoY, Ry2.sum.XtoY, p.all, n1, n2, test.details, alpha)
    
    if (H0.idx.XtoY) {
      type.I.err.idx.ind <- (R2.LM.res.XtoY.ind$pval < alpha)
      type.I.err.idx.sum <- (R2.LM.res.XtoY.sum$pval < alpha)
      R2.LM.res.XtoY <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (R2.LM.res.XtoY.ind$pval < alpha)
      power.idx.sum <- (R2.LM.res.XtoY.sum$pval < alpha)
      R2.LM.res.XtoY <- c(power.idx.ind, power.idx.sum)
    }
    
    # Y to X
    R2.LM.res.YtoX.ind <- R2Test(Ry2.YtoX, Rx2.ind.YtoX, p.all, n1, n2, test.details, alpha)
    R2.LM.res.YtoX.sum <- R2Test(Ry2.YtoX, Rx2.sum.YtoX, p.all, n1, n2, test.details, alpha)
    if (H0.idx.YtoX) {
      type.I.err.idx.ind <- (R2.LM.res.YtoX.ind$pval < alpha)
      type.I.err.idx.sum <- (R2.LM.res.YtoX.sum$pval < alpha)
      R2.LM.res.YtoX <- c(type.I.err.idx.ind, type.I.err.idx.sum)
    }else{
      power.idx.ind <- (R2.LM.res.YtoX.ind$pval < alpha)
      power.idx.sum <- (R2.LM.res.YtoX.sum$pval < alpha)
      R2.LM.res.YtoX <- c(power.idx.ind, power.idx.sum)
    }
  }else{
    R2.LM.res.XtoY <- NULL
    R2.LM.res.YtoX <- NULL
    R2.LM.res.XtoY.ind <- NULL
    R2.LM.res.XtoY.sum <- NULL
    R2.LM.res.YtoX.ind <- NULL
    R2.LM.res.YtoX.sum <- NULL
  }
  
  if ("R2-GLM" %in% test.type) {
    # X to Y 
    
    if (simuData.attrs$causal.type == "XtoY") {
      mdl2 <- glm(Y2 ~ Z2, family = binomial(link = "logit"))
      mdl2.null <- glm(Y2 ~ 1, family = binomial(link = "logit"))
      R2.glm.R <- 1 - logLik(mdl2)/logLik(mdl2.null)
      yfitted <- mdl2$fitted.values
      R2.glm.ss <- glmR2_func(Y2, yfitted)[1]
      ## logistic model
      test.details <- list(type = "GLM", y = Y2, yhat = yfitted)
      # X to Y
      R2.GLM.res.R <- R2Test(Rx2.XtoY, R2.glm.R, p.all, n1, n2, test.details, alpha)
      R2.GLM.res.ss <- R2Test(Rx2.XtoY, R2.glm.ss, p.all, n1, n2, test.details, alpha)
      
      if (H0.idx.XtoY) {
        type.I.err.idx.ind <- (R2.GLM.res.R$pval < alpha)
        type.I.err.idx.sum <- (R2.GLM.res.ss$pval < alpha)
        R2.GLM.res <- c(type.I.err.idx.ind, type.I.err.idx.sum)
      }else{
        power.idx.ind <- (R2.GLM.res.R$pval < alpha)
        power.idx.sum <- (R2.GLM.res.ss$pval < alpha)
        R2.GLM.res <- c(power.idx.ind, power.idx.sum)
      }
    }else if (simuData.attrs$causal.type == "YtoX"){
      # Y to X
      mdl2 <- glm(X2 ~ Z2, family = binomial(link = "logit"))
      mdl2.null <- glm(X2 ~ 1, family = binomial(link = "logit"))
      R2.glm.R <- 1 - logLik(mdl2)/logLik(mdl2.null)
      xfitted <- mdl2$fitted.values
      R2.glm.ss <- glmR2_func(X2, xfitted)[1]
      ## logistic model
      test.details <- list(type = "GLM", y = X2, yhat = xfitted)
      
      R2.GLM.res.R <- R2Test(Ry2.YtoX, R2.glm.R, p.all, n1, n2, test.details, alpha)
      R2.GLM.res.ss <- R2Test(Ry2.YtoX, R2.glm.ss, p.all, n1, n2, test.details, alpha)
      
      if (H0.idx.YtoX) {
        type.I.err.R <- (R2.GLM.res.R$pval < alpha)
        type.I.err.ss <- (R2.GLM.res.ss$pval < alpha)
        R2.GLM.res <- c(type.I.err.R, type.I.err.ss)
      }else{
        power.idx.R <- (R2.GLM.res.R$pval < alpha)
        power.idx.ss <- (R2.GLM.res.ss$pval < alpha)
        R2.GLM.res <- c(power.idx.R, power.idx.ss)
      }
    }
  }else{
    R2.glm.R <- NULL
    R2.glm.ss <- NULL
    R2.GLM.res <- NULL
    R2.GLM.res.R <- NULL
    R2.GLM.res.ss <- NULL
  }
  
  # Steiger
  if ("Steiger" %in% test.type) {
    ## X to Y
    rx <- cor(X1, Z1)
    ry <- cor(Y2, Z2)
    Steiger.XtoY <- Steiger(rx, ry, n1, n2, alpha)
    Steiger.res.XtoY <- (Steiger.XtoY$pval < alpha)
    
    ## Y to X
    ry <- cor(Y1, Z1)
    rx <- cor(X2, Z2)
    Steiger.YtoX <- Steiger(ry, rx, n1, n2, alpha)
    Steiger.res.YtoX <- (Steiger.YtoX$pval < alpha)
  }else{
    Steiger.res.XtoY <- NULL
    Steiger.res.YtoX <- NULL
    Steiger.XtoY <- NULL
    Steiger.YtoX <- NULL
  }
  
  if ("CDRatio" %in% test.type) {
    ## X to Y
    rhoxz <- cor(X1, Z1)
    rhoyz <- cor(Y2, Z2)
    CDratio.XY <- CDratio(Z1, rhoxz, rhoyz, n1, n2, alpha)
    kyx.CI <- CDratio.XY$CIs.YX
    kxy.CI <- CDratio.XY$CIs.XY
    term1 <- (kyx.CI[2] < 0 & kyx.CI[1] > -1) || (kyx.CI[2] > 0 & kyx.CI[1] < 1)
    term2 <- (kxy.CI[1] > 1 ||kxy.CI[2] < -1)
    CDratio.res.XtoY <- term1 & term2
    
    ## Y to X
    rhoyz.2 <- cor(Y1, Z1)
    rhoxz.2 <- cor(X2, Z2)
    CDratio.YX <- CDratio(Z1, rhoyz.2, rhoxz.2, n1, n2, alpha)
    kyx.CI <- CDratio.XY$CIs.YX
    kxy.CI <- CDratio.XY$CIs.XY
    term1 <- (kyx.CI[2] < 0 & kyx.CI[1] > -1) || (kyx.CI[2] > 0 & kyx.CI[1] < 1)
    term2 <- (kxy.CI[1] > 1 ||kxy.CI[2] < -1)
    CDratio.res.YtoX <- term1 & term2
  }
  
  if ("R" %in% test.type) {
    ## X to Y
    Rx <- 
      R.res.XtoY <- RTest(Rx, Ry, n1, n2, alpha)
    
    ## Y to X
    R.res.YtoX <- RTest(Ry.2, Rx.2, n1, n2, alpha)
  }
  
  test.res.XtoY <- c(F.res.XtoY, F2.res.XtoY, R2.LM.res.XtoY, Steiger.res.XtoY, R2.GLM.res, CDratio.res.XtoY)
  test.res.YtoX <- c(F.res.YtoX, F2.res.XtoY, R2.LM.res.YtoX, Steiger.res.YtoX, R2.GLM.res, CDratio.res.YtoX)
  if ("R2-GLM" %in% test.type) {
    test.res.names <- c(sapply(test.type[!(test.type %in% c("Steiger", "R2-GLM", "CDRatio"))], function(x) paste(x, c("ind", "sum"), sep = ".")))
    test.res.names <- c(test.res.names, "Steiger", "R2.GLM.R", "R2.GLM.ss", "CDRatio")
  }else{
    test.res.names <- c(sapply(test.type[!(test.type %in% c("Steiger", "CDRatio"))], function(x) paste(x, c("ind", "sum"), sep = ".")))
    test.res.names <- c(test.res.names, "Steiger", "CDRatio")
  }
  
  return(list(causal.type = causal.type, 
              H0.idx = c(H0.idx.XtoY, H0.idx.YtoX),
              test.XtoY = test.res.XtoY,
              test.YtoX = test.res.YtoX, 
              test.names = test.res.names,
              r.square = c(Rx2.XtoY, Ry2.ind.XtoY, Ry2.sum.XtoY, Ry2.YtoX, Rx2.ind.YtoX, Rx2.sum.YtoX, R2.glm.R, R2.glm.ss), 
              test.res.ls = list(F.XtoY.ind = F.res.XtoY.ind, 
                                 F.XtoY.sum = F.res.XtoY.sum, 
                                 F.YtoX.ind = F.res.YtoX.ind, 
                                 F.YtoX.sum = F.res.YtoX.sum, 
                                 F2.XtoY.ind = F2.res.XtoY.ind, 
                                 F2.XtoY.sum = F2.res.XtoY.sum, 
                                 F2.YtoX.ind = F2.res.YtoX.ind, 
                                 F2.YtoX.sum = F2.res.YtoX.sum, 
                                 R2.LM.XtoY.ind = R2.LM.res.XtoY.ind, 
                                 R2.LM.XtoY.sum = R2.LM.res.XtoY.sum, 
                                 R2.LM.YtoX.ind = R2.LM.res.YtoX.ind, 
                                 R2.LM.YtoX.sum = R2.LM.res.YtoX.sum, 
                                 R2.GLM.R = R2.GLM.res.R,
                                 R2.GLM.ss = R2.GLM.res.ss, 
                                 Steiger.XtoY = Steiger.XtoY, 
                                 Steiger.YtoX = Steiger.YtoX, 
                                 CDratio.XtoY = CDratio.XY, 
                                 CDratio.YtoX = CDratio.YX)))
}

simu_Data_Normal <- function(n1, r12, p.XYU, corZ, alpha.YtoX, alpha.ZtoX, alpha.UtoX, beta.XtoY, beta.ZtoY, beta.UtoY, eta.ZtoU, error.X = 1, error.Y = 1, error.U = sqrt(2), seed = 1025) {
  p  <- sum(p.XYU, na.rm = T)
  n2 <- n1 * r12
  n  <- n1 + n2
  
  set.seed(seed)
  # generate effect of Z to X
  if (!is.na(p.XYU[1])) {
    p.snp.X <- p.XYU[1]
    if (is.na(alpha.ZtoX)) {
      alpha.effect <- runif(p.snp.X, 0.2, 0.3) * (rbinom(p.snp.X, 1, 0.5) * 2 - 1)
    }else if(length(alpha.ZtoX) == p.snp.X) {
      alpha.effect <- alpha.ZtoX
    }else if(length(alpha.ZtoX) == 2) {
      alpha.effect <- runif(p.snp.X, alpha.ZtoX[1], alpha.ZtoX[2]) * (rbinom(p.snp.X, 1, 0.5) * 2 - 1)
    }
  }else{
    p.snp.X <- 0
    alpha.effect <- NULL
  }
  
  # generate effect of Z to Y
  if (!is.na(p.XYU[2])) {
    p.snp.Y <- p.XYU[2]
    if (is.na(beta.ZtoY)) {
      beta.effect <- runif(p.snp.Y, 0.2, 0.3) * (rbinom(p.snp.Y, 1, 0.5) * 2 - 1)
    }else if(length(beta.ZtoY) == p.snp.Y) {
      beta.effect <- alpha.ZtoX
    }else if(length(beta.ZtoY) == 2) {
      beta.effect <- runif(p.snp.Y, beta.ZtoY[1], beta.ZtoY[2]) * (rbinom(p.snp.Y, 1, 0.5) * 2 - 1)
    }
  }else{
    p.snp.Y <- 0
    beta.effect <- NULL
  }
  
  # generate effect of invalid Zs
  if (!is.na(p.XYU[3])) {
    p.snp.U <- p.XYU[3]
    if (is.na(alpha.UtoX)){
      gamma.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U,1,0.5) * 2 - 1)
    }else if(length(alpha.UtoX) == 2){
      gamma.effect <- runif(p.snp.U,alpha.UtoX[1],alpha.UtoX[2])*(rbinom(p.snp.U,1,0.5)*2-1)
    }else if(length(alpha.UtoX) == 1){
      gamma.effect <- alpha.UtoX*(rbinom(p.snp.U,1,0.5)*2-1)
    }else{
      gamma.effect <- alpha.UtoX
    }
    
    if (is.na(beta.UtoX)){
      eta.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else if(length(beta.UtoX) == 2){
      eta.effect <- runif(p.snp.U, beta.UtoX[1],beta.UtoX[2]) * (rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else if(length(beta.UtoX) == 1){
      eta.effect <- beta.UtoX*(rbinom(p.snp.U, 1, 0.5) * 2 - 1)
    }else{
      eta.effect <- beta.UtoX
    }
    
    if (is.na(eta.ZtoU)) {
      xi.effect <- runif(p.snp.U, -0.2, 0.2)
    }else{
      xi.effect <- runif(p.snp.U, eta.ZtoU[1], eta.ZtoU[2])
    }
    p.all <- p.snp.X + p.snp.Y + p.snp.U
    
    # Z <- matrix(rnorm(n * p.all), ncol = p.all, byrow = T)
    Z <- rmvnorm(n, sigma = corZ)
    U <- Z %*% matrix(c(rep(0,p.snp.X),rep(0,p.snp.Y),xi.effect), ncol = 1) + rnorm(n, 0, error.U)
  }else{
    p.snp.U <- 0
    U <- rep(0, n)
    p.all <- p.snp.X + p.snp.Y
    gamma.effect <- NULL 
    eta.effect <- NULL
    xi.effect <- NULL
    Z <- matrix(rnorm(n * p.all), ncol = p.all, byrow = T)
  }
  
  eps.X <- rnorm(n, 0, error.X)
  eps.Y <- rnorm(n, 0, error.Y)
  
  X <- 0
  Y <- 0
  
  if (!is.null(alpha.effect)) {
    X <- Z[, 1:p.snp.X] %*% matrix(alpha.effect, ncol = 1)
    Y <- Z[, 1:p.snp.X] %*% matrix(beta.XtoY * alpha.effect, ncol = 1)
  }
  
  if (!is.null(beta.effect)) {
    X <- X + Z[, (p.snp.X + 1):(p.snp.X + p.snp.Y)] %*% matrix(alpha.YtoX * beta.effect, ncol = 1)
    Y <- Y + Z[, (p.snp.X + 1):(p.snp.X + p.snp.Y)] %*% matrix(beta.effect, ncol = 1)
  }
  
  if (!is.null(gamma.effect)) {
    if(!is.null(eta.effect)){
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(gamma.effect + alpha.YtoX * eta.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(beta.XtoY * gamma.effect + eta.effect, ncol = 1)
    }else{
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(gamma.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(beta.XtoY * gamma.effect, ncol = 1)
    }
  }else{
    if(!is.null(eta.effect)){
      X <- X + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(alpha.YtoX * eta.effect, ncol = 1)
      Y <- Y + Z[, (p.snp.X + p.snp.Y + 1):p.all] %*% matrix(eta.effect, ncol = 1)
    }
  }
  
  X = (X + (1 + alpha.YtoX) * U + eps.X + alpha.YtoX * eps.Y)/(1 - beta.XtoY * alpha.YtoX)
  Y = (Y + (1 + beta.XtoY) * U + beta.XtoY * eps.X + eps.Y)/(1 - beta.XtoY * alpha.YtoX)
  
  
  Z1 <- scale(Z[1:n1, ], scale = F)
  X1 <- scale(as.matrix(X[1:n1]), scale = F)
  Y1 <- scale(as.matrix(Y[1:n1]), scale = F)
  
  Z2 <- scale(Z[(n1 + 1):n, ], scale = F)
  X2 <- scale(as.matrix(X[(n1 + 1):n]), scale = F)
  Y2 <- scale(as.matrix(Y[(n1 + 1):n]), scale = F)
  
  ## get summary statistics
  ### for the first half
  corZ1 <- cor(Z1, Z1)
  z1Tz1 <- t(Z1) %*% Z1
  b_X <- as.numeric(1/diag(z1Tz1) * (t(Z1) %*% X1))
  rep_X1 <- X1[, rep(1, p.all)]
  se_X <- ifelse(p == 1, sqrt(colSums((X1 - Z1 * b_X)^2)/(n-2)/diag(z1Tz1)), sqrt(colSums((rep_X1 - Z1 %*% diag(b_X))^2)/(n-2)/diag(z1Tz1)))
  
  b_Y = as.numeric(1/diag(z1Tz1) * (t(Z1) %*% Y1))
  rep_Y1 = Y1[,rep(1,p.all)]
  se_Y = ifelse(p == 1,  
                sqrt(colSums((Y1 - Z1 * b_Y)^2)/(n-2)/diag(z1Tz1)), 
                sqrt(colSums((rep_Y1 - Z1 %*% diag(b_Y))^2)/(n-2)/diag(z1Tz1))
  )
  
  ### for the second half
  corZ2 <- cor(Z2, Z2)
  z2Tz2 <- t(Z2) %*% Z2
  b_X2 <- as.numeric(1/diag(z2Tz2) * (t(Z2) %*% X2))
  rep_X2 <- X2[, rep(1, p.all)]
  se_X2 <- ifelse(p == 1, sqrt(colSums((X2 - Z2 * b_X2)^2)/(n-2)/diag(z2Tz2)), sqrt(colSums((rep_X2 - Z2 %*% diag(b_X2))^2)/(n-2)/diag(z2Tz2)))
  
  b_Y2 = as.numeric(1/diag(z2Tz2) * (t(Z2) %*% Y2))
  rep_Y2 = Y2[ ,rep(1,p.all)]
  se_Y2 = ifelse(p == 1, 
                 sqrt(colSums((Y2 - Z2 * b_Y)^2)/(n-2)/diag(z2Tz2)),
                 sqrt(colSums((rep_Y2 - Z2 %*% diag(b_Y))^2)/(n-2)/diag(z2Tz2)))
  
  ## determine whether it's null or alternative 
  ## when testing X -> Y: H0: Rx2 <= Ry2 H1: Rx2 > Ry2
  ## when testing Y -> X: H0: Ry2 <= Rx2 H1: Ry2 > Rx2  
  if (alpha.YtoX == 0) {
    if (beta.XtoY == 0) {
      causal.type = "NonCausal"
    }else{
      causal.type = "XtoY"
    }
  }else{
    if (beta.XtoY == 0) {
      causal.type = "YtoX"
    }else{
      causal.type = "BiDirect"
    }
  }
  
  return(list(stage1_df = list(Y = Y1, X = X1, Z = Z1), 
              stage2_df = list(Y = Y2, X = X2, Z = Z2), 
              stat1 = list(corZ = corZ1, zTz = z1Tz1, b_X = b_X, se_X = se_X, b_Y = b_Y, se_Y = se_Y), 
              stat2 = list(corZ = corZ2, zTz = z2Tz2, b_X = b_X2, se_X = se_X2, b_Y = b_Y2, se_Y = se_Y2), 
              causal.type = causal.type))
} 

unbiasR2 <- function(R2, n, p) {
  res <- max(0, min(1, 1 - (n - 3) / (n - p - 1) * (1 - R2) * ifelse(1-R2 < 0.5, hypergeo_analytical_taylor((n - p - 1) / 2, 1 - R2), hypergeo_analytical((n - p - 1) / 2, 1 - R2))))
  return(res)
}

loop_func_XtoY_ind_linear <- function(n1, r12, test.type, simuData.attrs, alpha = 0.05, MAF.mag = NULL, seed = 1025, cor_cutoff = 0.8, K = 30) {
  # get sample size
  n2 <- n1 * r12
  n <- n1 + n2
  test.details <- list(type = "LM")
  
  ## generate data
  simuData.type = simuData.attrs$type
  if (simuData.type == "Bi") {
    p.XYU <- simuData.attrs$p.XYU
    coef.multi.X <- simuData.attrs$coef.multi.X
    coef.multi.Y <- simuData.attrs$coef.multi.Y
    multi.error <- simuData.attrs$multi.error
    mdl2.type <- simuData.attrs$mdl2.type
    p.all <- sum(p.XYU, na.rm = T)
    
    # get coefficients
    alpha.YtoX <- coef.multi.X[[1]]
    alpha.ZtoX <- coef.multi.X[[2]]
    alpha.UtoX <- coef.multi.X[[3]]
    
    beta.XtoY <- coef.multi.Y[[1]]
    beta.ZtoY <- coef.multi.Y[[2]]
    beta.UtoY <- coef.multi.Y[[3]]
    
    error.X <- multi.error[1]
    error.Y <- multi.error[2]
    error.U <- multi.error[3]
    
    eta.ZtoU <- NA
    simuData <- simu_Data(n1, r12, p.XYU, alpha.YtoX, alpha.ZtoX, alpha.UtoX, beta.XtoY, beta.ZtoY, beta.UtoY, eta.ZtoU, error.X, error.Y, error.U, MAF.mag, seed)
  }else if(simuData.type == "Normal") {
    p.XYU <- simuData.attrs$p.XYU
    coef.multi.X <- simuData.attrs$coef.multi.X
    coef.multi.Y <- simuData.attrs$coef.multi.Y
    multi.error <- simuData.attrs$multi.error
    mdl2.type <- simuData.attrs$mdl2.type
    corZ <- simuData.attrs$corZ
    p.all <- sum(p.XYU, na.rm = T)
    
    # get coefficients
    alpha.YtoX <- coef.multi.X[[1]]
    alpha.ZtoX <- coef.multi.X[[2]]
    alpha.UtoX <- coef.multi.X[[3]]
    
    beta.XtoY <- coef.multi.Y[[1]]
    beta.ZtoY <- coef.multi.Y[[2]]
    beta.UtoY <- coef.multi.Y[[3]]
    
    error.X <- multi.error[1]
    error.Y <- multi.error[2]
    error.U <- multi.error[3]
    simuData <- simu_Data_Normal(n1, r12, p.XYU, corZ, alpha.YtoX, alpha.ZtoX, alpha.UtoX, beta.XtoY, beta.ZtoY, beta.UtoY, eta.ZtoU, error.X, error.Y, error.U, seed)
  }else if (simuData.type == "ADNI") {
    clearn_result <- simuData.attrs$clean_result
    effect.XYU <- simuData.attrs$effect.XYU
    sel_SNPs <- simuData.attrs$sel_SNPs
    null_Case <- simuData.attrs$null_Case
    k <- simuData.attrs$k
    causal.type <- simuData.attrs$causal.type
    simuData <- simu_ADNIdata(n1, r12, clean_result, effect.XYU, sel_SNPs, k, causal.type, null_Case, seed)
      
    p.all <- dim(simuData$stage1_df$Z)[2]
  }else if(simuData.type == "Real") {
    simuData <- simuData.attrs$Data
    p.all <- dim(simuData$stage1_df$Z)[2]
  }
  
  stage1_df <- simuData$stage1_df
  stage2_df <- simuData$stage2_df
  stat1 <- simuData$stat1
  stat2 <- simuData$stat2
  
  ## Get datasets
  X1 <- stage1_df$X
  Y1 <- stage1_df$Y
  Z1 <- stage1_df$Z
  
  X2 <- stage2_df$X
  Y2 <- stage2_df$Y
  Z2 <- stage2_df$Z
  
  # determine H0 and H1 
  causal.type <- simuData$causal.type
  H0.XtoY.type <- c("NonCausal", "YtoX")
  H1.XtoY.type <- c("XtoY", "BiDirect")
  H0.YtoX.type <- c("NonCausal", "XtoY")
  H1.YtoX.type <- c("YtoX", "BiDirect")
  
  H0.idx.XtoY <- (causal.type %in% H0.XtoY.type)
  H0.idx.YtoX <- (causal.type %in% H0.YtoX.type)
  
  ##### Stage I
  ## prune 
  cor_SNPs <- abs(cor(Z1))
  cor_SNPs <- (cor_SNPs < cor_cutoff)^2
  diag(cor_SNPs) <- 1
  i <- 1
  idx_sel <- NULL
  bY2 <- stat2$b_Y
  seY2 <- stat2$se_Y
  while (i < nrow(cor_SNPs)) {
    ind <- which(cor_SNPs[i, ] == 1)
    cor_SNPs <- as.matrix(cor_SNPs[ind, ind])
    Z1 <- Z1[, ind]
    Z2 <- Z2[, ind]
    bY2 <- bY2[ind]
    i <- i + 1 
  }
  
  if (ncol(cor_SNPs) > K) {
    idx <- order(abs(cor(X1, Z1)), decreasing = T)[1:K]
    Z1 <- scale(Z1[, idx])
    Z2 <- scale(Z2[, idx])
    bY2 <- bY2[idx]
  }
  p.after.sel <- dim(Z1)[2]
  
  ## linear model
  lm1.XtoY <- lm(X ~ ., data = data.frame(X = X1, Z = Z1))
  ## lm1.YtoX <- lm(Y ~ ., data = data.frame(Y = Y1, Z = Z1))
  
  Rx2.XtoY <- unbiasR2(summary(lm1.XtoY)$r.squared, n1, p.after.sel)
  Rx2.XtoY <- ifelse(Rx2.XtoY == 1, 1 - 1e-12, Rx2.XtoY)
  Rx2.XtoY <- ifelse(Rx2.XtoY == 0, 1e-12, Rx2.XtoY)
  ## Ry2.YtoX <- unbiasR2(summary(lm1.YtoX)$r.squared, n1, p.after.sel)
  Rx.XtoY <- sqrt(Rx2.XtoY)
  ## Ry.YtoX <- sqrt(Ry2.YtoX)

  rxz.XtoY <- cor(X1, Z1)
  ## rxz.YtoX <- cor(X2, Z2)
  ryz.XtoY <- cor(Y2, Z2)
  ## ryz.YtoX <- cor(Y1, Z1)
  
  ##### Stage II
  lm2.XtoY <- lm(Y ~ ., data = data.frame(Y = Y2, Z = Z2))
  ## lm2.YtoX <- lm(X ~ ., data = data.frame(X = X2, Z = Z2))
  
  Ry2.XtoY <- unbiasR2(summary(lm2.XtoY)$r.squared, n2, p.after.sel)
  Ry2.XtoY <- ifelse(Ry2.XtoY == 1, 1 - 1e-12, Ry2.XtoY)
  Ry2.XtoY <- ifelse(Ry2.XtoY == 0, 1e-12, Ry2.XtoY)
  ## Rx2.YtoX <- unbiasR2(summary(lm2.YtoX)$r.squared, n2, p.after.sel)
  Ry.XtoY <- sqrt(Ry2.XtoY)
  ## Rx.YtoX <- sqrt(Rx2.YtoX)
  

  ## test if R == 0 or R2 == 0
  test.rx2 <- pnorm(Rx2.XtoY / sqrt(2 * var_R2(Rx2.XtoY, n1, p.after.sel))) > 1 - alpha
  test.ry2 <- pnorm(Ry2.XtoY / sqrt(2 * var_R2(Ry2.XtoY, n2, p.after.sel))) > 1 - alpha
  # test.rx2 <- (pbeta(Rx2.XtoY, (p.after.sel - 1) / 2, (n1 - p.after.sel) / 2) > 1 - alpha)
  # test.ry2 <- (pbeta(Ry2.XtoY, (p.after.sel - 1) / 2, (n2 - p.after.sel) / 2) > 1 - alpha)
  # Frx <- (n1 - p.after.sel) * Rx2.XtoY / ((p.after.sel - 1) * (1 - Rx2.XtoY))
  # Fry <- (n2 - p.after.sel) * Ry2.XtoY / ((p.after.sel - 1) * (1 - Ry2.XtoY))
  
  Frx <- (n1 - p.after.sel) * summary(lm1.XtoY)$r.squared / ((p.after.sel - 1) * (1 - summary(lm1.XtoY)$r.squared))
  Fry <- (n2 - p.after.sel) * summary(lm2.XtoY)$r.squared / ((p.after.sel - 1) * (1 - summary(lm2.XtoY)$r.squared))
  test.rx <- Frx > qf(1 - alpha, p.after.sel - 1, n1 - p.after.sel)
  test.ry <- Fry > qf(1 - alpha, p.after.sel - 1, n2 - p.after.sel)
  
  
  # test.rx <- (1 - pnorm(Rx.XtoY / sqrt((1 - Rx2.XtoY)^2 / n1)) < alpha)
  # test.ry <- (1 - pnorm(Ry.XtoY / sqrt((1 - Ry2.XtoY)^2 / n2)) < alpha)
  # test.rx <- (Rx.XtoY - qnorm(1 - alpha /2) * sqrt((1 - Rx2.XtoY)^2 / n1) > 0)
  # test.ry <- (Ry.XtoY - qnorm(1 - alpha /2) * sqrt((1 - Ry2.XtoY)^2 / n2) > 0)
  
  ## F test
  F.res.XtoY <- fTest(Rx2.XtoY, Ry2.XtoY, p.after.sel, n1, n2, alpha = alpha)
  ## F.res.YtoX <- fTest(Ry2.YtoX, Rx2.YtoX, p.after.sel, n1, n2, alpha = alpha)
  
  # term1 <- (F.res.XtoY$CI[1] > 0)
  # term2 <- (F.res.YtoX$CI[1] < 0)
  # if(term1 & test.ry2) {
    # F.idx <- "XtoY"
  # }else{
    # term1 <- (F.res.XtoY$CI[1] < 0)
    # term2 <- (F.res.YtoX$CI[1] > 0)
    # if (term1 & test.ry2) {
      # F.idx <- "YtoX"
    # }else{
      # F.idx <- "NULL"
    # }
  # }
  
  if((F.res.XtoY$CI[1] > 0) & (test.ry) & (test.rx)) {
    F.XtoY.idx <- "XtoY"
  }else if((F.res.XtoY$CI[2] < 0) & (test.ry) & (test.rx)) {
    F.XtoY.idx <- "YtoX"
  }else{
    F.XtoY.idx <- "NULL"
  }
  
  # if (F.res.YtoX$CI[1] > 0) {
    # F.YtoX.idx <- "YtoX"
  # }else if (F.res.YtoX$CI[2] < 0) {
    # F.YtoX.idx <- "XtoY"
  # }else{
    # F.YtoX.idx <- "NULL"
  # }
  
  F.res <- c(F.res.XtoY$test.stat, F.res.XtoY$pval, F.res.XtoY$CI, F.XtoY.idx)
             # F.res.YtoX$test.stat, F.res.YtoX$pval, F.res.YtoX$CI, F.YtoX.idx, F.idx)

  ## R2 test 
  R2.res.XtoY <- R2Test(Rx2.XtoY, Ry2.XtoY, p.after.sel, n1, n2, test.details, alpha)
  # R2.res.YtoX <- R2Test(Ry2.YtoX, Rx2.YtoX, p.after.sel, n1, n2, test.details, alpha)
  
  # term1 <- (R2.res.XtoY$pval < alpha) & (R2.res.XtoY$CI[1] > 0)
  # term2 <- (R2.res.YtoX$pval > 1 - alpha) & (R2.res.YtoX$CI[2] < 0)
  # if (term1 & term2) {
    # R2.idx <- "XtoY"
  # }else{
    # term1 <- (R2.res.XtoY$pval > 1 - alpha) & (R2.res.XtoY$CI[2] < 0) 
    # term2 <- (R2.res.YtoX$pval < alpha) & (R2.res.YtoX$CI[1] > 0)
    # if (term1 & term2) {
      # R2.idx <- "YtoX"
    # }else{
      # R2.idx <- "NULL"
    # }
  # }
  
  if ((R2.res.XtoY$pval < alpha) & (R2.res.XtoY$CI[1] > 0) & (test.ry2) & (test.rx2)) {
    R2.XtoY.idx <- "XtoY"
  }else if((R2.res.XtoY$pval > 1 - alpha) & (R2.res.XtoY$CI[2] < 0) & (test.ry2) & (test.rx2)) {
    R2.XtoY.idx <- "YtoX"
  }else{
    R2.XtoY.idx <- "NULL"
  }

  # if ((R2.res.YtoX$pval < alpha) & (R2.res.YtoX$CI[1] > 0)) {
    # R2.YtoX.idx <- "YtoX"
  # }else if((R2.res.YtoX$pval > 1 - alpha) & (R2.res.YtoX$CI[2] < 0)) {
    # R2.YtoX.idx <- "XtoY"
  # }else{
    # R2.YtoX.idx <- "NULL"
  # }
  
  R2.res <- c(R2.res.XtoY$test.stat, R2.res.XtoY$pval, R2.res.XtoY$CI, R2.XtoY.idx)
              # R2.res.YtoX$test.stat, R2.res.YtoX$pval, R2.res.YtoX$CI, R2.YtoX.idx, R2.idx)
  
  ## R test 
  R.res.XtoY <- RTest(Rx.XtoY, Ry.XtoY, n1, n2, alpha)
  # R.res.YtoX <- RTest(Ry.YtoX, Rx.YtoX, n1, n2, alpha)
  
  ### get direction by testing results
  # term1 <- (R.res.XtoY$pval < alpha) & (R.res.XtoY$CI[1] > 0) 
  # term2 <- (R.res.YtoX$pval > 1 - alpha) & (R.res.YtoX$CI[2] < 0)
  # if (term1 & term2) {
    # R.idx <- "XtoY"
  # }else{
    # term1 <- (R.res.XtoY$pval > 1 - alpha) & (R.res.XtoY$CI[2] < 0) 
    # term2 <- (R.res.YtoX$pval < alpha) & (R.res.YtoX$CI[1] > 0)
    # if (term1 & term2) {
      # R.idx <- "YtoX"
    # }else{
      # R.idx <- "NULL"
    # }
  # }
  
  
  if ((R.res.XtoY$pval < alpha) & (R.res.XtoY$CI[1] > 0) & test.ry & test.rx){
    R.XtoY.idx <- "XtoY"
  }else if ((R.res.XtoY$pval > 1 - alpha) & (R.res.XtoY$CI[2] < 0) & test.ry & test.rx) {
    R.XtoY.idx <- "YtoX"
  }else{
    R.XtoY.idx <- "NULL"
  }

  # if ((R.res.YtoX$pval < alpha) & (R.res.YtoX$CI[1] > 0)){
    # R.YtoX.idx <- "YtoX"
  # }else if ((R.res.YtoX$pval > 1 - alpha) & (R.res.YtoX$CI[2] < 0)) {
    # R.YtoX.idx <- "XtoY"
  # }else{
    # R.YtoX.idx <- "NULL"
  # }
  
  R.res <- c(R.res.XtoY$test.stat, R.res.XtoY$pval, R.res.XtoY$CI, R.XtoY.idx)
             # R.res.YtoX$test.stat, R.res.YtoX$pval, R.res.YtoX$CI, R.YtoX.idx, R.idx)
  
  ## Steiger test
  ### X to Y
  Steiger.res.XtoY <- Steiger(rxz.XtoY, ryz.XtoY, n1, n2, alpha)
  
  if ((Steiger.res.XtoY$pval < alpha) & (Steiger.res.XtoY$test.stat > 0)) {
    Steiger.XtoY.idx <- "XtoY"
  }else if ((Steiger.res.XtoY$pval > 1 - alpha) & (Steiger.res.XtoY$test.stat < 0)) {
    Steiger.XtoY.idx <- "YtoX"
  }else{
    Steiger.XtoY.idx <- "NULL"
  }
  
  ### Y to X
  # Steiger.res.YtoX <- Steiger(ryz.YtoX, rxz.YtoX, n1, n2, alpha)
  
  # if ((Steiger.res.YtoX$pval < alpha) & (Steiger.res.YtoX$test.stat > 0)) {
    # Steiger.YtoX.idx <- "YtoX"
  # }else if ((Steiger.res.YtoX$pval > 1 - alpha) & (Steiger.res.YtoX$test.stat < 0)) {
    # Steiger.YtoX.idx <- "XtoY"
  # }else{
    # Steiger.YtoX.idx <- "NULL"
  # }
  
  Steiger.res <- c(Steiger.res.XtoY$test.stat, Steiger.res.XtoY$pval, Steiger.res.XtoY$CI, Steiger.XtoY.idx)
  #                  Steiger.res.YtoX$test.stat, Steiger.res.YtoX$pval, Steiger.res.YtoX$CI, Steiger.YtoX.idx)
  
  ## CD ratio test
  ### X to Y
  CDratio.res.XtoY <- CDratio(Z1, rxz.XtoY, ryz.XtoY, n1, n2, alpha)
  kyx.CI <- CDratio.res.XtoY$CIs.YX
  kxy.CI <- CDratio.res.XtoY$CIs.XY
  ### if it's X to Y: CIxy outside [-1, 1] and CIyx inside [-1, 0) or (0, 1]
  term1 <- (kyx.CI[2] < 0 & kyx.CI[1] > -1) || (kyx.CI[1] > 0 & kyx.CI[2] < 1)
  term2 <- (kxy.CI[1] > 1 ||kxy.CI[2] < -1)
  if (term1 & term2) {
    CDratio.XtoY.idx  <- "XtoY"
  }else{
    term1 <- (kxy.CI[2] < 0 & kxy.CI[1] > -1) || (kxy.CI[1] > 0 & kxy.CI[2] < 1)
    term2 <- (kyx.CI[1] > 1 ||kyx.CI[2] < -1)
    if (term1 & term2) {
      CDratio.XtoY.idx  <- "YtoX"
    }else{
      CDratio.XtoY.idx  <- "NULL"
    }
  }
  
  ### Y to X
  # CDratio.res.YtoX <- CDratio(Z1, rxz.YtoX, ryz.YtoX, n2, n1, alpha)
  # kyx.CI <- CDratio.res.YtoX$CIs.YX
  # kxy.CI <- CDratio.res.YtoX$CIs.XY
  ### if it's X to Y: CIxy outside [-1, 1] and CIyx inside [-1, 0) or (0, 1]
  # term1 <- (kyx.CI[2] < 0 & kyx.CI[1] > -1) || (kyx.CI[1] > 0 & kyx.CI[2] < 1)
  # term2 <- (kxy.CI[1] > 1 ||kxy.CI[2] < -1)
  # if (term1 & term2) {
    # CDratio.YtoX.idx  <- "XtoY"
  # }else{
    # term1 <- (kxy.CI[2] < 0 & kxy.CI[1] > -1) || (kxy.CI[1] > 0 & kxy.CI[2] < 1)
    # term2 <- (kyx.CI[1] > 1 ||kyx.CI[2] < -1)
    # if (term1 & term2) {
      # CDratio.YtoX.idx  <- "YtoX"
    # }else{
      # CDratio.YtoX.idx  <- "NULL"
    # }
  # }
  
  CDratio.res <- c(CDratio.res.XtoY$test.stat, CDratio.res.XtoY$CIs.YX, CDratio.res.XtoY$CIs.XY, CDratio.XtoY.idx)
                 # rbind(c(CDratio.res.XtoY$test.stat, CDratio.res.XtoY$CIs.YX, CDratio.res.XtoY$CIs.XY, CDratio.XtoY.idx), 
                 #       c(CDratio.res.YtoX$test.stat, CDratio.res.YtoX$CIs.YX, CDratio.res.YtoX$CIs.XY, CDratio.YtoX.idx))
  
  CDratio.res <- matrix(CDratio.res, nrow = 1)
  colnames(CDratio.res) <- c("K.YX", "K.XY", "CIs.YX.l", "CIs.YX.u", "CIs.XY.l", "CIs.XY.u", "Index")
  # rownames(CDratio.res) <- c("XtoY", "YtoX")
  

  test.arr <- rbind(F.res, R.res, R2.res, Steiger.res)
  colnames(test.arr) <- c("XtoY.test.stat", "XtoY.pval", "XtoY.CI.l", "XtoY.CI.u", "XtoY.idx") 
                          # "YtoX.test.stat", "YtoX.pval", "YtoX.CI.l", "YtoX.CI.u", "YtoX.idx", "idx")
  rownames(test.arr) <- c("F", "R", "R2", "Steiger")


  # return(list(causal.type = causal.type, 
  #             H0.idx = c(H0.idx.XtoY, H0.idx.YtoX),
  #             r.squared.x = c(Rx2.XtoY, Rx2.YtoX), 
  #             r.squared.y = c(Ry2.XtoY, Ry2.YtoX), 
  #             r.x = c(Rx.XtoY, Rx.YtoX),
  #             r.y = c(Ry.XtoY, Ry.YtoX), 
  #             test.res = test.arr, 
  #             CDratio.res = CDratio.res))
  
  return(list(causal.type = causal.type, 
              H0.idx = c(H0.idx.XtoY, H0.idx.YtoX),
              r.squared = c(Rx2.XtoY, Ry2.XtoY), 
              r = c(Rx.XtoY, Ry.XtoY),
              test.res = test.arr, 
              CDratio.res = CDratio.res))
}



simu_func_XtoY_ind_linear <- function(M, n1, r12, test.type, simuData.attrs, alpha = 0.05, HyperOR = "XtoY", MAF.mag = NULL, seed = 1025) {
  set.seed(seed)
  seeds_vec <- sample(100 * M, M)
  res_all <- pblapply(1:M, function(i) loop_func_XtoY_ind_linear(n1, r12, test.type, simuData.attrs, alpha, MAF.mag, seeds_vec[i]))
  
  causal.type <- res_all[[1]]$causal.type
  H0.idx <- res_all[[1]]$H0.idx
  
  
  RandR2.arr <- t(sapply(res_all, function(x) c(x$r, x$r.squared)))
  colnames(RandR2.arr) <- c("Rx.XtoY", "Ry.XtoY", "Rx2.XtoY", "Ry2.XtoY")
    
  # r.squared.arr <- t(sapply(res_all, function(x) c(x$r.squared.x, x$r.squared.y)))
  # colnames(r.squared.arr) <- c("Rx2.XtoY", "Rx2.YtoX", "Ry2.XtoY", "Ry2.YtoX")
  
  # r.arr <- t(sapply(res_all, function(x) c(x$r.x, x$r.y)))
  # colnames(r.arr) <- c("Rx.XtoY", "Rx.YtoX", "Ry.XtoY", "Ry.YtoX")
  
  ## Get power or type I error 
  cat("True causal relationship type = ", causal.type, "\n")
  
  index.arr <- t(sapply(1:M, function(i) c(res_all[[i]]$test.res[, 5], res_all[[i]]$CDratio.res[7])))
  colnames(index.arr) <- c("F.XtoY.idx", "R.XtoY.idx", "R2.XtoY.idx", "Steiger.XtoY.idx", "CDratio.XtoY.idx")
                           # "F.YtoX.idx", "R.YtoX.idx", "R2.YtoX.idx", "Steiger.YtoX.idx",
                           # "F.idx", "R.idx", "R2.idx", "CDratio.XtoY.idx", "CDratio.YtoX.idx")
  
  if(HyperOR == "XtoY") {
    cat("H0: Rx2 <= Ry2 v.s. H1: Rx2 > Ry2\n")
    if (H0.idx[1]){
      cat("Current setting is under H0. Type I error is going to be calculated.\n")
      power.idx <- "TypeIerror"
    }else{
      cat("Current setting is under H1. Power is going to be calculated.\n")
      power.idx <- "Power"
    }
  }else{
    cat("H0: Ry2 <= Rx2 v.s. H1: Ry2 > Rx2\n")
    if (H0.idx[2]){
      cat("Current setting is under H0. Type I error is going to be calculated.\n")
      power.idx <- "TypeIerror"
    }else{
      cat("Current setting is under H1. Power is going to be calculated.\n")
      power.idx <- "Power"
    }
  }
  
  prob.tab <- apply(index.arr, 2, function(x) c(mean(x == "XtoY"), mean(x == "YtoX"), mean(x == "NULL")))
  
  rownames(prob.tab) <- c("XtoY", "YtoX", "NULL")
  print(prob.tab)
  
  # get statistics and p-values
  
  all.test.res <- t(sapply(res_all, function(x) c(x$test.res, x$CDratio.res)))
  
  col.names.all <- c(t(sapply(c("F", "R", "R2", "Steiger"), function(x) paste(x, c("XtoY.test.stat", "XtoY.pval", "XtoY.CI.l", "XtoY.CI.u", "XtoY.idx"), sep = "."))))
                                                                                 # "YtoX.test.stat", "YtoX.pval", "YtoX.CI.l", "YtoX.CI.u", "YtoX.idx", "idx"), sep = "."))))
  col.names.all.2 <- c("CDratio.XtoY.K.YX", "CDratio.XtoY.K.XY", "CDratio.XtoY.CIs.YX.l", "CDratio.XtoY.CIs.YX.u", "CDratio.XtoY.CIs.XY.l", "CDratio.XtoY.CIs.XY.u", "CDratio.XtoY.idx")
                       # "CDratio.YtoX.K.YX", "CDratio.YtoX.K.XY", "CDratio.YtoX.CIs.YX.l", "CDratio.YtoX.CIs.YX.u", "CDratio.YtoX.CIs.XY.l", "CDratio.YtoX.CIs.XY.u", "CDratio.YtoX.idx")
  
  colnames(all.test.res) <- c(col.names.all, col.names.all.2)
  idx <- which(!grepl("idx", colnames(all.test.res)))
  all.test.res <- data.frame(all.test.res)
  for (idxi in idx) {
    all.test.res[, idxi] <- as.numeric(all.test.res[, idxi]) 
  }
  # all.test.res <- all.test.res[, which(colnames(all.test.res) != "Steiger.idx")]
  all.test.res$power.idx <- power.idx
  all.test.res$causal.type <- causal.type
  all.test.res$HyperOR <- HyperOR

  # idx <- which(!grepl("idx", colnames(all.test.res)))
  
  cat("## R and R.squared ##\n")
  print(summary(RandR2.arr))
  
  # cat("## R.squared ##\n")
  # print(summary(r.squared.arr))
  
  # cat("## R ##\n")
  # print(summary(r.arr))
  
  # cat("## test.statistics ##\n")
  # print(summary(all.test.res[, idx], na.rm = T))
  
  # return(list(TypeIorIIError = prob.tab, 
  #             test.res = all.test.res, 
  #             r.squared = r.squared.arr, 
  #             r = r.arr
  # ))
  return(list(TypeIorIIError = prob.tab, 
              test.res = all.test.res, 
              RandR2 = RandR2.arr
  ))
}



simu_func <- function(M, n1, r12, test.type, simuData.attrs, alpha = 0.05, HyperOR = "XtoY", MAF.mag = NULL, seed = 1025) {
  set.seed(seed)
  seeds_vec <- sample(100 * M, M)
  res_all <- pblapply(1:M, function(i) loop_func(n1, r12, test.type, simuData.attrs, alpha, MAF.mag, seeds_vec[i]))
  
  causal.type <- res_all[[1]]$causal.type
  H0.idx <- res_all[[1]]$H0.idx
  test.names <- res_all[[1]]$test.names
  
  ## Get power or type I error 
  cat("True causal relationship type = ", causal.type, "\n")
  if(HyperOR == "XtoY") {
    cat("H0: Rx2 <= Ry2 v.s. H1: Rx2 > Ry2\n")
    if (H0.idx[1]){
      cat("Current setting is under H0. Type I error is going to be calculated.\n")
    }else{
      cat("Current setting is under H1. Power is going to be calculated.\n")
    }
    
    test.res <- colMeans(t(sapply(1:M, function(i) res_all[[i]]$test.XtoY)), na.rm = T)
  }else{
    cat("H0: Ry2 <= Rx2 v.s. H1: Ry2 > Rx2\n")
    if (H0.idx[2]){
      cat("Current setting is under H0. Type I error is going to be calculated.\n")
    }else{
      cat("Current setting is under H1. Power is going to be calculated.\n")
    }
    
    test.res <- colMeans(t(sapply(1:M, function(i) res_all[[i]]$test.YtoX)), na.rm = T)
  }
  test.res <- matrix(test.res, nrow = 1)
  colnames(test.res) <- test.names
  print(test.res)
  # print(test.res[1, c(1,3,5,7)])
  
  # get statistics and p-values
  
  if ("R2-GLM" %in% test.type) {
    r.squared.arr <- array(NA, dim = c(M, 8), dimnames = list(1:M, c("Rx2.XtoY", "Ry2.ind.XtoY", "Ry2.sum.XtoY", "Ry2.YtoX", "Rx2.ind.YtoX", "Rx2.sum.YtoX", "R2.glm.R", "R2.glm.ss")))
    stats.arr <- array(NA, dim = c(M, 2 * length(test.names)), 
                       dimnames = list(1:M, c(paste(test.names[1:(length(test.names) - 2)], "XtoY", sep = "."), paste(test.names[1:(length(test.names) - 2)], "YtoX", sep = "."), "R2.GLM.R", "R2.GLM.ss", "CDRatio.XtoY.Kyx", "CDRatio.XtoY.Kxy", "CDRatio.YtoX.Kxy", "CDRatio.YtoX.Kyx")))
    pvals.arr <- array(NA, dim = c(M, 2 * (length(test.names) - 1) - 2), 
                       dimnames = list(1:M, c(paste(test.names[1:(length(test.names) - 2)], "XtoY", sep = "."), paste(test.names[1:(length(test.names) - 2)], "YtoX", sep = "."), "R2.GLM.R", "R2.GLM.ss")))
    
    for (i in 1:M) {
      stats.arr[i, ] <- c(res_all[[i]]$test.res.ls$F.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F2.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F2.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F2.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F2.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$Steiger.XtoY$test.stat, 
                          res_all[[i]]$test.res.ls$Steiger.YtoX$test.stat,
                          res_all[[i]]$test.res.ls$R2.GLM.R$test.stat, 
                          res_all[[i]]$test.res.ls$R2.GLM.ss$test.stat, 
                          res_all[[i]]$test.res.ls$CDratio.XtoY$test.stat, 
                          res_all[[i]]$test.res.ls$CDratio.YtoX$test.stat
      )
      pvals.arr[i, ] <- c(res_all[[i]]$test.res.ls$F.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$F.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$F.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$F.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$F2.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$F2.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$F2.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$F2.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$Steiger.XtoY$pval, 
                          res_all[[i]]$test.res.ls$Steiger.YtoX$pval,
                          res_all[[i]]$test.res.ls$R2.GLM.R$pval, 
                          res_all[[i]]$test.res.ls$R2.GLM.ss$pval)
      r.squared.arr[i, ] <- res_all[[i]]$r.square
    }
  }else{
    stats.arr <- array(NA, dim = c(M, 2 * length(test.names) + 2), 
                       dimnames = list(1:M, c(paste(test.names[1:(length(test.names)-1)], "XtoY", sep = "."), paste(test.names[1:(length(test.names)-1)], "YtoX", sep = "."), "CDRatio.XtoY.Kyx", "CDRatio.XtoY.Kxy", "CDRatio.YtoX.Kxy", "CDRatio.YtoX.Kyx")))
    pvals.arr <-  array(NA, dim = c(M, 2 * length(test.names) - 2), 
                        dimnames = list(1:M, c(paste(test.names[1:(length(test.names)-1)], "XtoY", sep = "."), paste(test.names[1:(length(test.names)-1)], "YtoX", sep = "."))))
    r.squared.arr <- array(NA, dim = c(M, 6), dimnames = list(1:M, c("Rx2.XtoY", "Ry2.ind.XtoY", "Ry2.sum.XtoY", "Ry2.YtoX", "Rx2.ind.YtoX", "Rx2.sum.YtoX")))
    
    for (i in 1:M) {
      stats.arr[i, ] <- c(res_all[[i]]$test.res.ls$F.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F2.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F2.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$F2.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$F2.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.ind$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.sum$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.ind$test.stat, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.sum$test.stat, 
                          res_all[[i]]$test.res.ls$Steiger.XtoY$test.stat, 
                          res_all[[i]]$test.res.ls$Steiger.YtoX$test.stat,
                          res_all[[i]]$test.res.ls$CDratio.XtoY$test.stat, 
                          res_all[[i]]$test.res.ls$CDratio.YtoX$test.stat
      )
      pvals.arr[i, ] <- c(res_all[[i]]$test.res.ls$F.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$F.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$F.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$F.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$F2.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$F2.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$F2.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$F2.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.ind$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.XtoY.sum$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.ind$pval, 
                          res_all[[i]]$test.res.ls$R2.LM.YtoX.sum$pval, 
                          res_all[[i]]$test.res.ls$Steiger.XtoY$pval, 
                          res_all[[i]]$test.res.ls$Steiger.YtoX$pval)
      r.squared.arr[i, ] <- res_all[[i]]$r.square
    }
  }
  
  return(list(TypeIorIIError = test.res, 
              test.stat = stats.arr, 
              pvals = pvals.arr, 
              r.squared = r.squared.arr))
}


loop_func_1 <- function(n, p, r, theta_XtoY, error.X = 1, error.Y = 1, alpha = alpha, seed  = 1025, MAF = 0.3) {
  n1 <- n
  n2 <- n * r
  n <- n1 + n2
  
  set.seed(seed)
  U <- rnorm(n, mean = 0, sd = error.U)
  Z <- matrix(rbinom(n * p, 2, MAF), nrow = n, ncol = p)
  beta.ZtoX <- runif(p, 0.2, 0.3) * (rbinom(p, 1, 0.5) * 2 - 1)
  
  err.X <- rnorm(n, sd = error.X)
  err.Y <- rnorm(n, sd = error.Y)
  
  X <- Z %*% matrix(beta.ZtoX, ncol = 1) + U * alpha.UtoX + err.X
  f.mat <- matrix(NA, n, 5)
  f.mat[, 1] <- X * theta_XtoY[1]
  f.mat[, 2] <- X^2 * theta_XtoY[2]
  f.mat[, 3] <- X^3 * theta_XtoY[3]
  f.mat[, 4] <- X^4 * theta_XtoY[4]
  f.mat[, 5] <- ifelse(X <= 0, 0, log(X)) * theta_XtoY[5]
  Y <- rowSums(f.mat) + beta.UtoY * U + err.Y
  
  if(all(theta_XtoY == 0)) {
    causal.type <- "NonCasual"
  }else{
    causal.type <- "XtoY"
  }
  
  idx <- sample(n, n)
  idx1 <- idx[1:n1]
  idx2 <- idx[(n1 + 1):n]
  
  X1 <- X[idx1]
  X2 <- X[idx2]
  Y1 <- Y[idx1]
  Y2 <- Y[idx2]
  Z1 <- Z[idx1, ]
  Z2 <- Z[idx2, ]
  
  ############## 
  lm1 <- lm(X1 ~ poly(Z1))
  lm2 <- lm(Y2 ~ poly(Z2))
  
  Rx2.XtoY <- unbiasR2(summary(lm1)$r.squared, n1, p)
  Ry2.XtoY <- unbiasR2(summary(lm2)$r.squared, n2, p)
  Rx2.XtoY <- ifelse(Rx2.XtoY == 0, 1e-12, Rx2.XtoY)
  Ry2.XtoY <- ifelse(Ry2.XtoY == 0, 1e-12, Ry2.XtoY)
  
  test.rx2 <- pnorm(Rx2.XtoY / sqrt(var_R2(Rx2.XtoY, n1, p))) > 1 - alpha
  test.ry2 <- pnorm(Ry2.XtoY / sqrt(var_R2(Ry2.XtoY, n2, p))) > 1 - alpha
  
  Frx <- (n1 - p) * Rx2.XtoY / ((p - 1) * (1 - Rx2.XtoY))
  Fry <- (n2 - p) * Ry2.XtoY / ((p - 1) * (1 - Ry2.XtoY))
  test.rx <- Frx > qf(1 - alpha, p - 1, n1 - p)
  test.ry <- Fry > qf(1 - alpha, p - 1, n2 - p)
  
  ## F test
  F.res.XtoY <- fTest(Rx2.XtoY, Ry2.XtoY, p, n1, n2, alpha = alpha)
  
  if((F.res.XtoY$CI[1] > 0) & (test.ry) & (test.rx)) {
    F.XtoY.idx <- "XtoY"
  }else if((F.res.XtoY$CI[2] < 0) & (test.ry) & (test.rx)) {
    F.XtoY.idx <- "YtoX"
  }else{
    F.XtoY.idx <- "NULL"
  }
  
  F.res <- c(F.res.XtoY$test.stat, F.res.XtoY$pval, F.res.XtoY$CI, F.XtoY.idx)
  
  ## R2 test 
  test.details <- list(type = "LM")
  R2.res.XtoY <- R2Test(Rx2.XtoY, Ry2.XtoY, p, n1, n2, test.details, alpha)
  
  if ((R2.res.XtoY$pval < alpha) & (R2.res.XtoY$CI[1] > 0) & (test.ry2) & (test.rx2)) {
    R2.XtoY.idx <- "XtoY"
  }else if((R2.res.XtoY$pval > 1 - alpha) & (R2.res.XtoY$CI[2] < 0) & (test.ry2) & (test.rx2)) {
    R2.XtoY.idx <- "YtoX"
  }else{
    R2.XtoY.idx <- "NULL"
  }
  
  R2.res <- c(R2.res.XtoY$test.stat, R2.res.XtoY$pval, R2.res.XtoY$CI, R2.XtoY.idx)
  
  ## R test 
  R.res.XtoY <- RTest(sqrt(Rx2.XtoY), sqrt(Ry2.XtoY), n1, n2, alpha)
  
  if ((R.res.XtoY$pval < alpha) & (R.res.XtoY$CI[1] > 0) & test.ry & test.rx){
    R.XtoY.idx <- "XtoY"
  }else if ((R.res.XtoY$pval > 1 - alpha) & (R.res.XtoY$CI[2] < 0) & test.ry & test.rx) {
    R.XtoY.idx <- "YtoX"
  }else{
    R.XtoY.idx <- "NULL"
  }
  
  
  R.res <- c(R.res.XtoY$test.stat, R.res.XtoY$pval, R.res.XtoY$CI, R.XtoY.idx)
  
  ## Steiger test
  ### X to Y
  rxz.XtoY <- cor(X1, Z1)
  ryz.XtoY <- cor(Y2, Z2)
  Steiger.res.XtoY <- Steiger(rxz.XtoY, ryz.XtoY, n1, n2, alpha)
  
  if ((Steiger.res.XtoY$pval < alpha) & (Steiger.res.XtoY$test.stat > 0)) {
    Steiger.XtoY.idx <- "XtoY"
  }else if ((Steiger.res.XtoY$pval > 1 - alpha) & (Steiger.res.XtoY$test.stat < 0)) {
    Steiger.XtoY.idx <- "YtoX"
  }else{
    Steiger.XtoY.idx <- "NULL"
  }
  
  Steiger.res <- c(Steiger.res.XtoY$test.stat, Steiger.res.XtoY$pval, Steiger.res.XtoY$CI, Steiger.XtoY.idx)
  
  ## CD ratio test
  ### X to Y
  CDratio.res.XtoY <- CDratio(Z1, rxz.XtoY, ryz.XtoY, n1, n2, alpha)
  kyx.CI <- CDratio.res.XtoY$CIs.YX
  kxy.CI <- CDratio.res.XtoY$CIs.XY
  ### if it's X to Y: CIxy outside [-1, 1] and CIyx inside [-1, 0) or (0, 1]
  term1 <- (kyx.CI[2] < 0 & kyx.CI[1] > -1) || (kyx.CI[1] > 0 & kyx.CI[2] < 1)
  term2 <- (kxy.CI[1] > 1 ||kxy.CI[2] < -1)
  if (term1 & term2) {
    CDratio.XtoY.idx  <- "XtoY"
  }else{
    term1 <- (kxy.CI[2] < 0 & kxy.CI[1] > -1) || (kxy.CI[1] > 0 & kxy.CI[2] < 1)
    term2 <- (kyx.CI[1] > 1 ||kyx.CI[2] < -1)
    if (term1 & term2) {
      CDratio.XtoY.idx  <- "YtoX"
    }else{
      CDratio.XtoY.idx  <- "NULL"
    }
  }
  
  CDratio.res <- c(CDratio.res.XtoY$test.stat, CDratio.res.XtoY$CIs.YX, CDratio.res.XtoY$CIs.XY, CDratio.XtoY.idx)
  
  CDratio.res <- matrix(CDratio.res, nrow = 1)
  colnames(CDratio.res) <- c("K.YX", "K.XY", "CIs.YX.l", "CIs.YX.u", "CIs.XY.l", "CIs.XY.u", "Index")
  
  
  test.arr <- rbind(F.res, R.res, R2.res, Steiger.res)
  colnames(test.arr) <- c("XtoY.test.stat", "XtoY.pval", "XtoY.CI.l", "XtoY.CI.u", "XtoY.idx") 
  rownames(test.arr) <- c("F", "R", "R2", "Steiger")
  
  
  return(list(causal.type = causal.type, 
              r.squared = c(Rx2.XtoY, Ry2.XtoY, sqrt(var_R2(Rx2.XtoY, n1, p)), sqrt(var_R2(Ry2.XtoY, n2, p))), 
              test.res = test.arr, 
              CDratio.res = CDratio.res))
} 



simu_func_1 <- function(M, n, p, r, theta_XtoY, alpha = 0.05, seed = 1025) {
  set.seed(seed)
  seeds_vec <- sample(100 * M, M)
  res_all <- pblapply(1:M, function(i) loop_func_1(n, p, r, theta_XtoY, alpha = alpha, seed = seeds_vec[i]))
  
  causal.type <- res_all[[1]]$causal.type
  
  RandR2.arr <- t(sapply(res_all, function(x) x$r.squared))
  colnames(RandR2.arr) <- c("Rx2", "Ry2", "se.Rx2", "se.Ry2")
  RandR2.arr.vec <- c(colMeans(RandR2.arr[, 1:2]), apply(RandR2.arr[, 1:2], 2, sd), colMeans(RandR2.arr[, 3:4]))
  RandR2.arr.vec <- matrix(RandR2.arr.vec, nrow = 1)
  colnames(RandR2.arr.vec) <- c("ave.Rx2", "ave.Ry2", "sd.Rx2", "sd.Ry2", "se.Rx2", "se.Ry2")
  
  cat("True causal relationship type = ", causal.type, "\n")
  cat("n1 = ", n, " | n2 = ", n * r, " | p = ", p, " | Current theta_XtoY = ", theta_XtoY, "\n")
  
  index.arr <- t(sapply(1:M, function(i) c(res_all[[i]]$test.res[, 5], res_all[[i]]$CDratio.res[7])))
  colnames(index.arr) <- c("F.XtoY.idx", "R.XtoY.idx", "R2.XtoY.idx", "Steiger.XtoY.idx", "CDratio.XtoY.idx")
  
  prob.tab <- apply(index.arr, 2, function(x) c(mean(x == "XtoY"), mean(x == "YtoX"), mean(x == "NULL")))
  
  rownames(prob.tab) <- c("XtoY", "YtoX", "NULL")
  print(prob.tab)
  
  # get statistics and p-values
  
  all.test.res <- t(sapply(res_all, function(x) c(x$test.res, x$CDratio.res)))
  col.names.all <- c(t(sapply(c("F", "R", "R2", "Steiger"), function(x) paste(x, c("XtoY.test.stat", "XtoY.pval", "XtoY.CI.l", "XtoY.CI.u", "XtoY.idx"), sep = "."))))
  col.names.all.2 <- c("CDratio.XtoY.K.YX", "CDratio.XtoY.K.XY", "CDratio.XtoY.CIs.YX.l", "CDratio.XtoY.CIs.YX.u", "CDratio.XtoY.CIs.XY.l", "CDratio.XtoY.CIs.XY.u", "CDratio.XtoY.idx")
  
  colnames(all.test.res) <- c(col.names.all, col.names.all.2)
  idx <- which(!grepl("idx", colnames(all.test.res)))
  all.test.res <- data.frame(all.test.res)
  for (idxi in idx) {
    all.test.res[, idxi] <- as.numeric(all.test.res[, idxi]) 
  }
  
  all.test.res$causal.type <- causal.type
  
  cat("## R and R.squared ##\n")
  print(summary(RandR2.arr))
  
  cat("## R summary ##\n")
  print(RandR2.arr.vec)
  
  return(list(TypeIorIIError = prob.tab, 
              test.res = all.test.res, 
              RandR2 = RandR2.arr,
              summary.R = RandR2.arr.vec
  ))
}


