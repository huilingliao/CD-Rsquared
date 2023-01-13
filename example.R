################################### Data generation ###################################
n <- 50000 # sample size
p.snp.X <- 15 # number of SNPs in g_X
p.snp.Y <- 10 # number of SNPs in g_Y
p.snp.U <- 10 # number of SNPs in g_B
xi.range <- 0.2 # generate non-zero xi's from uniform distribution, implying correlated pleiotropy
theta.XtoY <- 0.1 # causal effect from X to Y
theta.YtoX <- 0   # causal effect from Y to X
seed <- 1025

p.all <- p.snp.X + p.snp.Y + p.snp.U
set.seed(seed)

# generate alpha
alpha.effect <- runif(p.snp.X, 0.2, 0.3) * (rbinom(p.snp.X, 1, 0.5)*2 - 1)

# generate beta
beta.effect <- runif(p.snp.Y, 0.2, 0.3) * (rbinom(p.snp.Y, 1, 0.5)*2 - 1)

# generate effects of invalid IVs
gamma.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U, 1, 0.5)*2 - 1)
eta.effect <- runif(p.snp.U, 0.2, 0.3) * (rbinom(p.snp.U, 1, 0.5)*2 - 1)
xi.effect <- runif(p.snp.U, -xi.range, xi.range)


# generate individual data
MAF <- rep(0.3, p.all)
Z <- matrix(rbinom(2 * n * p.all, 2, MAF), ncol = p.all, byrow = T)
U <- Z %*% c(rep(0, p.snp.X), rep(0, p.snp.Y), xi.effect) + rnorm(2 * n, 0, sqrt(2))
error_X <- rnorm(2 * n, 0, 1)
error_Y <- rnorm(2 * n, 0, 1)

X <- (Z %*% c(alpha.effect, theta.YtoX * beta.effect, gamma.effect + theta.YtoX * eta.effect) + (1 + theta.YtoX) * U + error_X + theta.YtoX * error_Y) / (1 - theta.XtoY * theta.YtoX)
Y <- (Z %*% c(theta.XtoY * alpha.effect, beta.effect, gamma.effect * theta.XtoY + eta.effect) + (1 + theta.XtoY) * U + error_Y + theta.XtoY * error_X) / (1 - theta.XtoY * theta.YtoX)

# generate two independent samples (Z1, X1, Y1) and (Z2, X2, Y2)
Z1 <- Z[1:n, ]
X1 <- as.matrix(X[1:n])
Y1 <- as.matrix(Y[1:n])

Z2 <- Z[(1 + n):(2 * n), ]
X2 <- as.matrix(X[(1 + n):(2 * n)])
Y2 <- as.matrix(Y[(1 + n):(2 * n)])

Z1 <- scale(Z1, scale = F)
X1 <- scale(X1, scale = F)
Y1 <- scale(Y1, scale = F)

Z2 <- scale(Z2, scale = F)
X2 <- scale(X2, scale = F)
Y2 <- scale(Y2, scale = F)

# get summary statistics for X and Y from marginal linear regression
z1Tz1 <- t(Z1) %*% Z1
b_X <- as.numeric(1 / diag(z1Tz1) * (t(Z1) %*% X1))
rep_X1 <- X1[, rep(1, p.all)]
se_X <- sqrt(colSums((rep_X1 - Z1 %*% diag(b_X))^2) / (n - 2)/diag(z1Tz1))

z2Tz2 <- t(Z2) %*% Z2
b_Y <- as.numeric(1 / diag(z2Tz2) * (t(Z2) %*% Y2))
rep_Y2 <- Y2[, rep(1, p.all)]
se_Y <- sqrt(colSums((rep_Y2 - Z2 %*% diag(b_Y))^2) / (n - 2)/diag(z2Tz2))

################################### Simulation ###################################

lm1 <- lm(X1 ~ Z1)
rx2 <- summary(lm1)$r.squared

lm2 <- lm(Y2 ~ Z2)
ry2 <- summary(lm2)$r.squared

type <- "LM"
yhat <- predict(lm2)

R2Test(rx2, ry2, p.all, n, n, alpha, type, Y2, yhat)