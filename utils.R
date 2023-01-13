# library(Rcpp)
# sourceCpp("D:/Projects/Casual-Rsquared/Simulation code/code/matrix_multi.cpp")

## data generation
#'@param oriSNP the small dataset for data generation
#'@param u.sd standard deviation for u
#'@param num.sample the total number of samples if sample.type == 'sample', the number of replicates if sample.type == 'duplicate'
#'@param sample.type choose from 'sample' (sample with replacement) or 'duplicate' (duplicate the whole dataset)
#'@param stage1_mdl fitted model for stage 1
#'@param stage2_mdl fitted model for stage 2
#'@param effect.size.U effect size for U
#'@param effect.size.X effect size for X
#'@param split.r sample split ratio, n2 = r x n1
genData <- function(oriSNP, u.sd, num.sample, sample.type = 'sample', mdl2.type, sigma_Y, 
                    stage1_mdl, stage2_mdl, stage2_mdl_null, effect.size.U, effect.size.X, split.r = 1, direction = "XtoY", setseed = 1025) {
  # get original dataset size
  n0 <- nrow(oriSNP)
  
  # sample U ~ N(0, hat_u_sd)
  set.seed(setseed)
  
  # sample with replacement to get SNPs
  if(sample.type == 'sample') {
    nall <- num.sample
    SNPs_idx <- sample(n0, num.sample, replace = TRUE)
    SNPs <- oriSNP[SNPs_idx, ]
  }else{
    SNPs <- NULL
    for (i in 1:num.sample) {
      SNPs <- rbind(SNPs, oriSNP)
    }
    nall <- num.sample * n0
  }
  U <- rnorm(num.all, sd = u.sd)
  
  mdl2 <- strsplit(mdl2.type, "-")[[1]][1]
  # generate X
  # include.idx <- sapply(strsplit(rownames(summary(stage1_mdl)$coefficients), "SNP")[-1], function(x) as.numeric(x[2]))
  # X <- as.matrix(cbind(1, SNPs[, include.idx])) %*% matrix(stage1_mdl$coefficients, ncol = 1) + U
  X <- predict(stage1_mdl, newdata = data.frame(SNPs)) + U
  
  # generate y
  if (mdl2 == "GLM") {
    if (effect.size.X == 0) {
      logit.p <- cbind(1, U) %*% matrix(stage2_mdl_null$coefficients * c(1, effect.size.U), ncol = 1) + rnorm(num.sample, 0, sigma_Y)
    }else{
      logit.p <- cbind(1, X, U) %*% matrix(stage2_mdl$coefficients * c(1, effect.size.X, effect.size.U), ncol = 1) + rnorm(num.sample, 0, sigma_Y)
    }
    p <- exp(logit.p) / (1 + exp(logit.p))
    y <- rbinom(num.sample, 1, p)
  }else if (mdl2 == "LM") {
    if (effect.size.X == 0) {
      y <- cbind(1, U) %*% matrix(stage2_mdl_null$coefficients * c(1, effect.size.U), ncol = 1) + rnorm(num.sample, 0, sigma_Y)
    }else{
      y <- cbind(1, X, U) %*% matrix(stage2_mdl$coefficients * c(1, effect.size.X, effect.size.U), ncol = 1) + rnorm(num.sample, 0, sigma_Y)
    }
  }

  # get datasets for two stages
  n1 <- floor(nall / (1 + split.r))
  idx1 <- sample(nall, n1)
  idx2 <- setdiff(1:nall, idx1)

  y1 <- scale(y[idx1], scale = FALSE)
  y2 <- scale(y[idx2], scale = FALSE)
  X1 <- scale(X[idx1], scale = FALSE)
  X2 <- scale(X[idx2], scale = FALSE)
  SNP1 <- scale(SNPs[idx1, ], scale = FALSE)
  SNP2 <- scale(SNPs[idx2, ], scale = FALSE)

  if (direction == "XtoY") {
    stage1.df <- data.frame(Y = y1, X = X1, SNP = SNP1)
    stage2.df <- data.frame(Y = y2, X = X2, SNP = SNP2)
  }else{
    stage1.df <- data.frame(Y = X1, X = y1, SNP = SNP1)
    stage2.df <- data.frame(Y = X2, X = y2, SNP = SNP2)
  }
  colnames(stage1.df)[3:ncol(stage1.df)] <- colnames(oriSNP)
  colnames(stage2.df)[3:ncol(stage2.df)] <- colnames(oriSNP)

  return(list(stage1.df = stage1.df,
              stage2.df = stage2.df))
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
fTest <- function(rx2, ry2, M, n1, n2, lbd1, lbd2, num.sample = 10000) {
  statx <- rx2 * (n1 - M - 1) / ((1 - rx2) * M)
  staty <- ry2 * (n2 - M - 1) / ((1 - ry2) * M)
  # nu1 <- max(0, (statx * M - M) / n1)
  # nu2 <- max(0, (staty * M - M) / n2)
  nu <- 0.5 * (lbd1 / n1 + lbd2 / n2)
  lbdx <- nu * n1 # nu1 * n1
  lbdy <- nu * n2 # nu2 * n2
  
  Tx <- rf(num.sample, M, n1 - M - 1, lbdx)
  Ty <- rf(num.sample, M, n2 - M - 1, lbdy)

  TnuX <- (Tx * M - M) / n1
  TnuY <- (Ty * M - M) / n2
  
  pval_xy <- mean(Tx < Ty)
  # pval_xy <- mean(statx - staty < Tx - Ty)
  # pval_yx <- mean(staty - statx < Ty - Tx)
  return(list(test.stat = statx - staty, pval = pval_xy))
}


## Testing depends on beta
#'@param rx2 scalar, the R-squared for stage 1 regression
#'@param ry2 scalar, the R-squared for stage 2 regression
#'@param M the number of SNPs
#'@param n1 the sample size for stage 1
#'@param n2 the sample size for stage 2
#'@param lbd1 noncentrality parameter for stage 1 regression
#'@param lbd2 noncentrality parameter for stage 2 regression
#'@param num.sample number of samples for bootstrap
#'

beta_mu_var <- function(alp, bt, lbd){
  uu <- c(alp + bt, alp + 1)
  ll <- c(alp, alp + bt + 1)
  z <- lbd/2
  # mu <- exp(-lbd/2) * gamma(alp + 1) / gamma(alp) * gamma(alp + bt) / gamma(alp + bt + 1) * genhypergeo(uu, ll, z)
  mu <- exp(-lbd/2) * alp / (alp + bt) * genhypergeo(uu, ll, z)
  uu <- c(alp + bt, alp + 2)
  ll <- c(alp, alp + bt + 2)
  sigma <- exp(-lbd/2) * alp * (alp + 1) / ((alp + bt) * (alp + bt + 1)) * genhypergeo(uu, ll, z) - mu^2
  if(exp(-lbd/2) == 0) {
    mu <- 0
    sigma <- 0
  }
  return(c(mu, sigma))
}


betaTest <- function(rx2, ry2, M, n1, n2, lbd1, lbd2, num.sample) {
  # nu <- 0.5 * (lbd1 / n1 + lbd2 / n2)
  # lbdx <- nu * n1
  # lbdy <- nu * n2
  
  # Tx <- rbeta(num.sample, M/2, (n1 - M - 1)/2, lbd1)
  # Ty <- rbeta(num.sample, M/2, (n2 - M - 1)/2, lbd2)
  
  # pval_xy <- mean(Tx < Ty)
  
  # pval_xy <- mean(rx2 - ry2 < Tx - Ty)
  # pval_yx <- mean(ry2 - rx2 < Ty - Tx)
  
  # nu1 <- rx2 /(1 - rx2)
  # nu2 <- ry2 /(1 - ry2)
  # lbdx <- nu1 * n1
  # lbdy <- nu2 * n2
  # Tx <- rbeta(num.sample, M/2, (n1 - M - 1)/2, lbdx)
  # Ty <- rbeta(num.sample, M/2, (n2 - M - 1)/2, lbdy)
  # pval_xy <- mean(Tx < Ty)
  # use normal
  var1 <- beta_mu_var(M/2, (n1 - M - 1)/2, lbd1)[2]
  var2 <- beta_mu_var(M/2, (n2 - M - 1)/2, lbd2)[2]
  zstat <- (rx2 - ry2) / sqrt(var1 + var2)
  pval_xy <- 1 - pnorm(zstat)
  return(list(test.stat = zstat, pval = pval_xy))
}

## Testing based on CLT for R2
#'@param y binary response
#'@param phat fitted logistic values
glmR2_func <- function(y, phat) {
  SST <- sum((y - mean(y))^2)
  SSE <- sum((y - phat)^2)
  Rss2 <- 1 - SSE / SST
  Rr2 <- sum((y - mean(y)) * (phat - mean(y)))^2 / (sum((y - mean(y))^2) * sum((phat - mean(y))^2))
  return(c(Rss2, Rr2))
}

#'@param rx2 scalar, the R-squared for stage 1 regression
#'@param ry2 scalar, the R-squared for stage 2 regression
#'@param M the number of SNPs
#'@param n1 the sample size for stage 1
#'@param n2 the sample size for stage 2
#'@param alpha significant level
#'@param type type of regression model
#'@param yhat fitted values
R2Test <- function(rx2, ry2, M, n1, n2, alpha, type, y, yhat) {
  if(type == "LM") {
    # var1 <- 4 * rx2 * (1 - rx2)^2 * (n1 - M - 1)^2 / ((n1^2 - 1) * (n1 + 3))
    # var2 <- 4 * ry2 * (1 - ry2)^2 * (n2 - M - 1)^2 / ((n2^2 - 1) * (n2 + 3))
    var1 <- 4 / n1 * rx2 * (1 - rx2)^2 * (1 -  (2 * M + 5) / n1)
    var2 <- 4 / n2 * ry2 * (1 - ry2)^2 * (1 -  (2 * M + 5) / n2)
    estsd <- sqrt(var1 + var2)
  }else{
    pp <- mean(yhat)
    V1 <- mean(yhat * (1 - yhat))
    V2 <- var(y)
    # V2 <- pp * (1 - pp)
    Z <- cbind(y, y * yhat, yhat^2)
    Sig <- cov(Z)
    c1 <- matrix(c((V2 - V1 + 2 * pp * V1) / V2^2, -2 / V2, 1 / V2), nrow = 1)
    estsd <- sqrt((c1 %*% Sig %*% t(c1))[1])
  }
  zstat <- (rx2 - ry2) / estsd
  pval_xy <- 1 - pnorm(zstat)

  # CIs <- c(rx2 - ry2 + qnorm(alpha/2) * estsd, rx2 - ry2 + qnorm(1- alpha/2) * estsd)
  return(list(test.stat = zstat, pval = pval_xy))
}


## Steiger's test
Steiger <- function(df1, df2) {
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  p <- ncol(df1) 
  
  # select the SNP for test statistic construction
  rx <- cor(df1$X, df1[, 3:p])
  ry <- cor(df2$Y, df2[, 3:p])
  rxy <- abs(rx) + abs(ry)
  idx <- which.max(rxy)
  
  # construct the test statistics
  zx <- 0.5 * log((1 + abs(rx[idx]))/(1 - abs(rx[idx])))
  zy <- 0.5 * log((1 + abs(ry[idx]))/(1 - abs(ry[idx])))
  z <- (zx - zy) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  pval_xy <- 1 - pnorm(z)
  
  return(list(test.stat = z, pval = pval_xy))
}


## RatioCD
Pfunc <- function(x, SigMat) {
  M <- length(x)
  Pmat <- matrix(NA, M + 1, M + 1)
  Pmat[1:M, 1:M] <- SigMat
  Pmat[M + 1, 1:M] <- x
  Pmat[1:M, M + 1] <- x
  Pmat[1 + M, 1 + M] <- 1
  return(Pmat)
}

generate_V <- function(SNP,rho)
{
  SNP = t(SNP)
  n = nrow(SNP)
  Sigma = rho[1:n,1:n]
  inv_Sigma = solve(Sigma,tol = 0)
  rho_X = as.matrix(rho[1:n,(n+1)])
  alpha =  inv_Sigma %*% rho_X
  e2 = 1 - t(rho_X) %*% inv_Sigma %*% rho_X
  
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

# ratioCD <- function(df1, df2, SigMat = NULL) {
#   n1 <- nrow(df1)
#   n2 <- nrow(df2)
#   p <- ncol(df1) 
#   
#   SNP1 <- df1[, 3:p]
#   SNP2 <- df2[, 3:p]
#   for (i in 1:(p-2)) {
#     SNP1[, i] <- scale(SNP1[, i])
#     SNP2[, i] <- scale(SNP2[, i])
#   }
#   X1 <- scale(df1$X)
#   Y2 <- scale(df2$Y)
#   # select the SNP for test statistic construction
#   rx <- cor(X1, SNP1)
#   ry <- cor(Y2, SNP2)
#   
#   if(is.null(SigMat)) {
#     SNPs <- as.matrix(rbind(SNP1, SNP2))
#     SigMat <- t(SNPs) %*% SNPs / (n1 + n2)
#   }
#   PxMat <- Pfunc(rx, SigMat)
#   PyMat <- pfunc(ry, SigMat)
#   
#   newX <- cbind(df1[, 3:p], df1$X)
#   newY <- cbind(df2[, 3:p], df2$Y)
#   
#   Vx <- generate_V(newX, PxMat)
#   Vy <- generate_V(newY, PyMat)
#   Jmat <- cbind(diag(1/rx), -diag(ry/rx^2))
#   combinedV <- rbind(cbind(Vy, matrix(0, ncol = p, nrow = p)))
# }


## one testing 
simu1_func <- function(oriSNP, uhat.sd, num.sample, sample.type = 'sample', 
                       stage1_mdl, stage2_mdl, stage2_mdl_null, sigma_Y, 
                       effect.size.U, effect.size.X, split.r, setseed = 1025,
                       mdl1.type = "LM",  mdl2.type = "LM-RI", test.type, num.pseudo = 10000) {
  # data generation
  df_lst <- genData(oriSNP, uhat.sd, num.sample, sample.type, mdl2.type, sigma_Y, 
                    stage1_mdl, stage2_mdl, stage2_mdl_null, effect.size.U, effect.size.X, split.r, setseed)
  df1 <- df_lst[[1]]
  df2 <- df_lst[[2]]
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  
  # model fitting
  ## 1st stage
  if (mdl1.type == "LM") {
    lm1.formula <- paste0("X ~ ", paste0(colnames(df1)[grepl("SNP", colnames(df1))], collapse = " + "))
    lm1.formula <- as.formula(lm1.formula)
    mdl1 <- lm(lm1.formula, data = df1)
  }else{
    print("Method not implemented. ")
  }
  
  SNPs <- df1[, 3:(ncol(df1))]
  M <- ncol(df1) - 2
  # Rx2 <- matrix(cor(df1$X, SNPs), nrow = 1) %*% solve(cor(SNPs, SNPs), tol = 1e-40) %*% matrix(cor(df1$X, SNPs), ncol = 1)
  Rx2 <- summary(mdl1)$r.squared
  lbdx <- n1 * Rx2 / (1 - Rx2)
  
  ## 2nd stage
  mdl2.type <- strsplit(mdl2.type, "-")[[1]]
  SNPs2 <- df2[, 3:(ncol(df2))]
  if (mdl2.type[2] == "RI") {
    uhat <- df2$Y - predict(mdl1, newdata = df2)
    if (mdl2.type[1] == "LM") {
      mdl2 <- lm(df2$Y ~ df2$X + uhat)
      Ry2 <- as.numeric(matrix(cor(df2$Y, SNPs2), nrow = 1) %*% solve(cor(SNPs2, SNPs2), tol = 1e-30) %*% matrix(cor(df2$Y, SNPs2), ncol = 1))
      # Ry2 <- summary(mdl2)$r.squared
      yhat <- mdl2$fitted.values - mdl2$coefficients[3] * uhat
      newY <- df2$Y - mdl2$coefficients[3] * uhat
    }else{
      mdl2 <- glm(df2$Y ~ df2$X + uhat, family = binomial(link = "logit"))
      # mdl2.2 <- glm(df2$Y ~ uhat, family = binomial(link = "logit"))
      # mdl2 <- glm(df2$Y ~ df2$X, family = binomial(link = "logit"))
      # yhat <- mdl2$fitted.values
      yhat <- 1 / (1 + exp( - mdl2$coefficients[1] - mdl2$coefficients[2] * df2$X))
      Ry2 <- glmR2_func(df2$Y, yhat)[2]
      # newY <- ?
      # Ry2 <- summary(mdl2.1)$r.squared - summary(mdl2.2)$r.squared
    }
  }else if (mdl2.type[2] == "PS") {
    xhat <- predict(mdl1, newdata = df2)
    if (mdl2.type[1] == "LM") {
      mdl2 <- lm(df2$Y ~ xhat)
      # Ry2 <- as.numeric(matrix(cor(df2$Y, SNPs2), nrow = 1) %*% solve(cor(SNPs2, SNPs2), tol = 1e-30) %*% matrix(cor(df2$Y, SNPs2), ncol = 1))
      Ry2 <- summary(mdl2)$r.squared
      yhat <- mdl2$fitted.values
      newY <- df2$Y
    }else{
      mdl2 <- glm(df2$Y ~ xhat, family = binomial(link = "logit"))
      yhat <- mdl2$fitted.values
      Ry2 <- glmR2_func(df2$Y, yhat)[2]
      newY <- df2$Y
    }
  }
  lbdy <- n2 * Ry2 / (1 - Ry2)
  
  # hypothesis testing
  if(length(test.type) == 1) {
    test_res <- switch (test.type,
                        "F" = fTest(Rx2, Ry2, M, n1, n2, lbdx, lbdy, num.pseudo), 
                        "Beta" = betaTest(Rx2, Ry2, M, n1, n2, lbdx, lbdy, num.pseudo), 
                        "R2" = R2Test(Rx2, Ry2, M, n1, n2, alpha, mdl2.type[1], newY, yhat),
                        "Steiger" = Steiger(df1, df2)
    )
    test.res <- c(test_res$test.stat, test_res$pval)
    res.cols <- c(paste0(test.type, "test.stat"), paste0(test.type, "pval_xy"))
  }else{
    test_res <- rep(NA, 2 * length(test.type))
    res.cols <- c(sapply(test.type, function(x) paste(x, c("test.stat", "pval_xy"))))
    for (k in 1:length(test.type)) {
      test_resk <- switch (test.type[k],
                           "F" = fTest(Rx2, Ry2, M, n1, n2, lbdx, lbdy, num.pseudo), 
                           "Beta" = betaTest(Rx2, Ry2, M, n1, n2, lbdx, lbdy, num.pseudo),
                           "R2" = R2Test(Rx2, Ry2, M, n1, n2, alpha, mdl2.type[1], df2$Y, yhat), 
                           "Steiger" = Steiger(df1, df2)
                           )
      test_res[(k - 1) * 2 + 1] = test_resk$test.stat
      test_res[(k - 1) * 2 + 2] = test_resk$pval
    }
    
  }
  
  return(list(r.squared = c(Rx2, Ry2), lambda = c(lbdx, lbdy), testRes = test_res, res.cols = res.cols))
}



## the complete simulation function

simu_func <- function(M, oriSNP, uhat.sd, num.sample, sample.type = 'sample', 
                      stage1_mdl, stage2_mdl, stage2_mdl_null, sigma_Y, 
                      effect.size.U, effect.size.X, split.r, setseed = 1025, 
                      mdl1.type = "lm", mdl2.type = "LM-RI", test.type = "R2", num.pseudo = 10000) {
  
  res <- pblapply(1:M, function(i) simu1_func(oriSNP, uhat.sd, num.sample, sample.type, 
                                            stage1_mdl, stage2_mdl, stage2_mdl_null, sigma_Y, 
                                            effect.size.U, effect.size.X, split.r, 
                                            setseed = setseed * i + 1992,
                                            mdl1.type,  mdl2.type, test.type, num.pseudo))
  r.squared <- t(sapply(res, function(x) x$r.squared)) # M x 2
  lambda.xy <- t(sapply(res, function(x) x$lambda)) # M x 2
  test.res <- t(sapply(res, function(x) x$testRes))
  colnames(test.res) <- res[[1]]$res.cols
  return(list(r.squared = r.squared, 
              lambda.xy = lambda.xy, 
              test.res = test.res))
}