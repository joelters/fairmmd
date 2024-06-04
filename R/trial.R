# ### Data ####
# rm(list = ls())
# # load_all()
# library(MASS)
# library(dplyr)
# n <- 5000
# k <- 10
# w <- 2
# X <- mvrnorm(n,rep(0,k),diag(k))
# b <- (-1)^c(1:k)
# noise <- 0.1
# # X1 <- rnorm(n)
# # X2 <- rnorm(n)
# # X3 <- rnorm(n)
# A <- runif(n) > 0.5
# # X <- data.frame(X1 = X1, X2 = X2, X3 = X3, S = A)
# eps <- rnorm(n,noise)
# Yl <- w*A + X%*%b + eps
# Y <- as.numeric(Yl > 0)
# FVs <- pnorm(w*A + X%*%b, mean = 0, sd = noise)
# X <- data.frame(X, S = A)
# df <- as_tibble(cbind(Y = Y,X, A = A))
# res <- NULL
# res[[1]] <- df
# res[[2]] <- FVs
# res[[3]] <- fair_infeasible(FVs = FVs, A = A, pac1 = A == 1, pac0 = A == 0)
# names(res) <- c("df", "FVt", "mmdt")
# res$mmdt
#
# PI <- fairmmd(Y,X,A,est_method = "Plugin", ML = "RF", fairness = "SP")
# D1 <- fairmmd(Y,X,A,est_method = "Debiased", CFit = FALSE, ML = "RF", fairness = "SP",
#               sterr = TRUE)
# D2 <- fairmmd(Y,X,A,est_method = "Debiased", CFit = TRUE, ML = "RF", fairness = "SP",
#               sterr = TRUE, verbose = TRUE)
# # DEB <- fairmmd(Y,X,A,est_method = "Debiased", ML = "RF", fairness = "SP")
# res$mmdt
# PI$FairMMD
# D1
# D2$fairMMD
