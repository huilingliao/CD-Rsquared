#####################  Load packages  #####################
source("D:/Projects/Casual-Rsquared/Simulation code/code/utils.R")
library(hypergeo)
library(pbapply)
pboptions(type = "txt", char = "=")

#####################  Load sample dataset  #####################

load("D:/Projects/Casual-Rsquared/twas_methods-master/sample_gene_expression_and_snp.Rdata")

colnames(cleaning_result)[4:33] = paste("SNP",1:30,sep="")
colnames(cleaning_result)[2] <- "y"
colnames(cleaning_result)[3] <- "X"

lm_formula = paste("X ~ ", paste("SNP", 1:30, sep = "",collapse =" + "))
lm_formula = as.formula(lm_formula)
lm1 = step(lm(lm_formula,data = cleaning_result), direction = "backward", trace = -1)
# include.idx <- sapply(strsplit(rownames(summary(lm1)$coefficients), "SNP")[-1], function(x) as.numeric(x[2]))

# Predict y and u
X_hat = predict(lm1)
u_hat = cleaning_result$X - X_hat
uhat.sd = sd(u_hat)

# # Linear regression for stage 2
lm2.1 <- lm(cleaning_result$y ~ u_hat)
lm2.2 <- lm(cleaning_result$y ~ cleaning_result$X + u_hat)

# Logistic Regression Only u_hat
glm2.1 = glm(cleaning_result$y ~ u_hat, family = binomial(link = "logit"))

# Logistic Regression with y and u_hat
glm2.2 = glm(cleaning_result$y ~ cleaning_result$X + u_hat,
             family = binomial(link = "logit"))


###################################################################

# num.replicates <- 1
# effect.size.X<- 1
# effect.size.U <- 1
# split.r <- 1
# M <- 10 # replicates
# stage1_mdl <- lm1
# mdl1.type <- "LM"
# mdl2.type <- "LM-PS"
# stage2_mdl <- lm2.2
# stage2_mdl_null <- lm2.1
# oriSNP <- cleaning_result[, 4:33]
# sample.type = 'sample'
# num.pseudo <- 100
# setseed <- 1025
# sigma.Y <- 0.1
# test.type <- c("F", "Beta", "R2", "Steiger")
# num.sample <- num.replicates * nrow(oriSNP)
# 
# simu_func(M, oriSNP, uhat.sd, num.sample, sample.type,
#           stage1_mdl, stage2_mdl, stage2_mdl_null, sigma.Y, 
#           effect.size.U, effect.size.X, split.r, setseed,
#           mdl1.type, mdl2.type, test.type, num.pseudo)


#####################  Parameter Settings  #####################


num.replicates <- c(1, 3, 5, 10, 15, 20)
effect.size.X<- c(0, 1, 2, 3, 5, 10)
effect.size.U <- c(1, 10, 20, 30, 40, 50)
split.r <- c(1, 2, 3, 5, 10)
sigma.Y <- c(0, 0.05, 0.1, 0.2, 0.5, 1, 2)
M <- 1000 # replicates
stage1_mdl <- lm1 
mdl1.type <- "LM"
# mdl2.type <- c("LM-RI", "LM-PS", "GLM-RI", "GLM-PS")
mdl2.type <- c("LM-RI", "LM-PS", "GLM-PS")
stage2_mdl <- list(lm2.2, lm2.2, glm2.2)
stage2_mdl_null <- list(lm2.1, lm2.1, glm2.1)
# stage2_mdl <- glm2.2
# stage2_mdl_null <- glm2.1
oriSNP <- cleaning_result[, 4:33]
sample.type = 'sample'
num.pseudo <- 10000
setseed <- 1025
test.type <- list(c("F", "Beta", "R2", "Steiger"), 
                  c("F", "Beta", "R2", "Steiger"), 
                  c("R2"))

#####################  Simulation  #####################

# sample size
res <- vector("list", length(num.replicates))
for (i in 1:length(num.replicates)) {
  # split.r
  res[[i]] <- vector("list", length(split.r))
  for (j in 1:length(split.r)) {
    # effect.size.X
    res[[i]][[j]] <- vector("list", length(effect.size.X))
    for (k in 1:length(effect.size.X)) {
      # effect.size.U
      res[[i]][[j]][[k]] <- vector("list", length(effect.size.U))
      for (m in 1:length(effect.size.U)) {
        # mdl2.type
        res[[i]][[j]][[k]][[m]] <-  vector("list", length(mdl2.type))
        for (p in 1:length(mdl2.type)) {
          # sigma.Y size
          res[[i]][[j]][[k]][[m]][[p]] <- vector("list", length(sigma.Y))
          for (qq in 1:length(sigma.Y)) {
            res[[i]][[j]][[k]][[m]][[p]][[qq]] <- simu_func(M, oriSNP, uhat.sd, num.replicates[i] * nrow(oriSNP), sample.type, 
                                                      stage1_mdl, stage2_mdl[[p]], stage2_mdl_null[[p]], sigma.Y[qq], 
                                                      effect.size.U[m], effect.size.X[k], split.r[j], setseed, 
                                                      mdl1.type, mdl2.type[p], test.type[[p]], num.pseudo)
          }
        }
      }
    }
  }
}
