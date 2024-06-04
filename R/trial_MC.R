# ### Data ####
# rm(list = ls())
# # load_all()
# library(MASS)
# library(dplyr)
# n <- 100
# k <- 10
# w <- 2
# b <- (-1)^c(1:k)
# noise <- 0.1
# R <- 2
#
# ntr <- 20000
# X <- mvrnorm(ntr,rep(0,k),diag(k))
# A <- runif(ntr) > 0.5
# eps <- rnorm(ntr,noise)
# Yl <- w*A + X%*%b + eps
# Y <- as.numeric(Yl > 0)
# FVs <- pnorm(w*A + X%*%b, mean = 0, sd = noise)
# thtr <- fair_infeasible(FVs = FVs, A = A, pac1 = A == 1, pac0 = A == 0)
#
# res <- lapply(1:R, function(u){
#   print(u)
#   X <- mvrnorm(n,rep(0,k),diag(k))
#   A <- runif(n) > 0.5
#   eps <- rnorm(n,noise)
#   Yl <- w*A + X%*%b + eps
#   Y <- as.numeric(Yl > 0)
#   FVs <- pnorm(w*A + X%*%b, mean = 0, sd = noise)
#   #IMPORTANT THAT a IS NOT CALLED A IN X
#   X <- data.frame(X, S = A)
#   PI <- fairmmd(Y,X,A,est_method = "Plugin", ML = "RF", fairness = "SP", sterr = FALSE)
#   D2 <- fairmmd(Y,X,A,est_method = "Debiased", CFit = TRUE, ML = "RF", fairness = "SP",
#                 sterr = TRUE)
#   se <- D2$fairMMD[2]
#   PI_rej <- abs(PI - thtr)/se >= qnorm(0.975)
#   D2_rej <- abs(D2 - thtr)/se >= qnorm(0.975)
#   data.frame(PI = PI[1], D2 = D2$fairMMD[1], PI_rej = PI_rej, D2_rej = D2_rej)
# })
# res2 <- do.call(rbind,res)
# colMeans(res2)
