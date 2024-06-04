mlfair <- function(X,
                   Y,
                   A,
                   pac1,
                   pac0,
                   relative = TRUE,
                   ML = c("Lasso", "Ridge", "RF", "CIF", "XGB","Logit_lasso"),
                   fitted_values = TRUE,
                   rf.ntree = 500,
                   rf.depth = NULL,
                   polynomial = polynomial,
                   weights = NULL){
  ML = match.arg(ML)
  #Estimate FVs
  m <- ML::MLest(X,
                 Y,
                 ML,
                 rf.cf.ntree = rf.ntree,
                 rf.depth = rf.depth,
                 FVs = TRUE)
  FVs <- m$FVs
  if(sum(FVs < 0 | FVs > 1) != 0){
    warning(paste(sum(FVs <= 0),"There are estimated probabilities outside unit interval"))
  }
  #Estimate plug in fairMMD

  n <- length(Y)
  n1 <- n - 1
  Ppac0 <- mean(pac0)
  Ppac1 <- mean(pac1)

  fairpi <- sapply(1:n1, function(u){
    j <- u + 1
    a <- (2/(n*n1))*sum((exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u])*(pac0[j:n]))/(Ppac0*Ppac0) +
               (exp(-0.5*(FVs[u] - FVs[j:n])^2)*pac1[u]*pac1[j:n])/(Ppac1*Ppac1) -
               (exp(-0.5*(FVs[u] - FVs[j:n])^2)*pac1[u]*pac0[j:n])/(Ppac1*Ppac0) -
               (exp(-0.5*(FVs[u] - FVs[j:n])^2)*pac0[u]*pac1[j:n])/(Ppac0*Ppac1))

  })
  fairpi <- sum(fairpi)
  res <- fairpi
  names(res) <- "Fair"

  if(relative == TRUE){
    D <- 2 - 2*exp(-0.5)
    fairpirel <- fairpi/D
    res <- as.matrix(c(fairpi, fairpirel))
    rownames(res) <- c("Fair","Fair_rel")
  }

  if (fitted_values == TRUE){
    return(list(fairMMD = res, FVs = FVs))
  }
  else{return(res)}
}

fair_deb <- function(Y, FVs, pac1, pac0, pacx1, pacx0, A){
  n <- length(Y)
  n1 <- n - 1
  Ppac0 <- mean(pac0)
  Ppac1 <- mean(pac1)
  b1 <- lapply(1:n1, function(u){
    j <- u + 1
    th00 <- (2/(n*n1))*sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u]*pac0[j:n] -
                   (FVs[u] - FVs[j:n])*pacx0[u]*pacx0[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n])))/(Ppac0*Ppac0)
    th11 <- (2/(n*n1))*sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u]*pac1[j:n] -
                   (FVs[u] - FVs[j:n])*pacx1[u]*pacx1[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n])))/(Ppac1*Ppac1)
    th10 <- (2/(n*n1))*sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u]*pac0[j:n] -
                   (FVs[u] - FVs[j:n])*pacx1[u]*pacx0[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n])))/(Ppac1*Ppac0)
    th01 <- (2/(n*n1))*sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u]*pac1[j:n] -
                   (FVs[u] - FVs[j:n])*pacx0[u]*pacx1[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n])))/(Ppac0*Ppac1)
    return(data.frame(th00 = th00, th11 = th11, th10 = th10, th01 = th01))
  })
  b1 <- do.call(rbind,b1)
  b1 <- colSums(b1)
  return(b1)
}

fairdnumsq <- function(Y1, Y2, FVs1,FVs2,A1,A2,
                       pac1_1, pac1_0, pacx1_1, pacx1_0,
                       pac2_1, pac2_0, pacx2_1, pacx2_0,
                       Ppac0,Ppac1){
  n1 <- length(FVs1)
  b1 <- lapply(1:n1, function(u){
    th00 <- sum(exp(-0.5*(FVs1[u] - FVs2)^2)*(pac1_0[u]*pac2_0 -
                   (FVs1[u] - FVs2)*pacx1_0[u]*pacx2_0*(Y1[u] - Y2 - FVs1[u] + FVs2))/(Ppac0*Ppac0))
    th11 <- sum(exp(-0.5*(FVs1[u] - FVs2)^2)*(pac1_1[u]*pac2_1 -
                   (FVs1[u] - FVs2)*pacx1_1[u]*pacx2_1*(Y1[u] - Y2 - FVs1[u] + FVs2))/(Ppac1*Ppac1))
    th10 <- sum(exp(-0.5*(FVs1[u] - FVs2)^2)*(pac1_1[u]*pac2_0 -
                   (FVs1[u] - FVs2)*pacx1_1[u]*pacx2_0*(Y1[u] - Y2 - FVs1[u] + FVs2))/(Ppac1*Ppac0))
    th01 <- sum(exp(-0.5*(FVs1[u] - FVs2)^2)*(pac1_0[u]*pac2_1 -
                   (FVs1[u] - FVs2)*pacx1_0[u]*pacx2_1*(Y1[u] - Y2 - FVs1[u] + FVs2))/(Ppac0*Ppac1))
    return(data.frame(th00 = th00, th11 = th11, th10 = th10, th01 = th01))
  })
  b1 <- do.call(rbind,b1)
  b1 <- colSums(b1)
  return(b1)
}

fairdnumtr <- function(Y, FVs, A,
                        pac1, pac0, pacx1, pacx0,
                        Ppac0, Ppac1){
  n <- length(FVs)
  n1 <- n-1
  b1 <- lapply(1:n1, function(u){
    j <- u + 1
    th00 <- sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u]*pac0[j:n] -
                   (FVs[u] - FVs[j:n])*pacx0[u]*pacx0[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n]))/(Ppac0*Ppac0))
    th11 <- sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u]*pac1[j:n] -
                   (FVs[u] - FVs[j:n])*pacx1[u]*pacx1[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n]))/(Ppac1*Ppac1))
    th10 <- sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u]*pac0[j:n] -
                   (FVs[u] - FVs[j:n])*pacx1[u]*pacx0[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n]))/(Ppac1*Ppac0))
    th01 <- sum(exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u]*pac1[j:n] -
                   (FVs[u] - FVs[j:n])*pacx0[u]*pacx1[j:n]*
                   (Y[u] - Y[j:n] - FVs[u] + FVs[j:n]))/(Ppac0*Ppac1))
    return(data.frame(th00 = th00, th11 = th11, th10 = th10, th01 = th01))
  })
  b1 <- do.call(rbind,b1)
  b1 <- colSums(b1)
  return(b1)
}


se_deb <- function(Y, FVs, A, pac1, pac0, pacx1, pacx0, th11, th00, th10, th01){
  n <- length(FVs)
  na0 <- sum(A == 0)
  na1 <- sum(A == 1)
  Sab <- rep(0,4)
  for (ii in 1:4){
    a <- as.numeric((ii == 1) + (ii == 2))
    b <- as.numeric((ii == 1) + (ii == 3))
    th <- th11*(ii == 1) + th00*(ii == 4) +
      th10*(ii == 2) + th01*(ii == 3)
    na = if (ii == 1 | ii == 2) na1 else na0
    nb = if (ii == 1 | ii == 3) na1 else na0
    paca = if (ii == 1 | ii == 2) pac1 else pac0
    pacb = if (ii == 1 | ii == 3) pac1 else pac0
    pacxa = if (ii == 1 | ii == 2) pacx1 else pacx0
    pacxb = if (ii == 1 | ii == 3) pac1 else pac0
    # browser()
    S <- sapply(1:n, function(u){
      aux <- sum((1/(n-1))*(exp(-0.5*(FVs[u] - FVs[-u])^2)*((paca[u]*pacb[-u] -
                              (FVs[u] - FVs[-u])*pacxa[u]*pacxb[-u]*(Y[u] - Y[-u] - FVs[u] + FVs[-u])))
                               - th*-(paca[u])*(pacb[-u])))
    })
    if (ii == 2){
      S10 <- S
    }
    else if (ii == 3){
      S01 <- S
    }
    else if (ii == 1){
      S11 <- S
    }
    else {
      S00 <- S
    }
  }
  # browser()
  SIG <- matrix(rep(0,16),4,4)
  for (i in 1:4){
    for (j in 1:4){
      SIG[i,j] <- (1/n)*sum(S11^2*(i == 1 & j == 1) +
                              S00^2*(i == 2 & j == 2) +
                              S10^2*(i == 3 & j == 3) +
                              S01^2*(i == 4 & j == 4) +
                              S11*S00*((i == 1 & j == 2) | i == 2 & j == 1) +
                              S11*S10*((i == 1 & j == 3) | i == 3 & j == 1) +
                              S11*S01*((i == 1 & j == 4) | i == 4 & j == 1) +
                              S00*S10*((i == 2 & j == 3) | i == 3 & j == 2) +
                              S00*S01*((i == 2 & j == 4) | i == 4 & j == 2) +
                              S10*S01*((i == 3 & j == 4) | i == 4 & j == 3))
    }
  }
  selc <- as.matrix(c(1,1,-1,-1))
  S <- t(selc) %*% SIG %*% selc
  B <- ((mean(paca))*(mean(pacb)))
  V = (4*S)/(B^2)
  se = sqrt(V/n)
  return(se)
}

SP <- function(df, npart){
  nn <- nrow(df)
  p <- base::split(sample(nn,nn,replace = FALSE),as.factor(1:npart))
  # p <- hyperSMURF::do.random.partition(nn, npart, seed = 0)
  dfsp <- NULL
  for (i in 1:npart){
    dfsp[[i]] <- dplyr::as_tibble(df[p[[i]],])
  }
  return(list(dfsp = dfsp, indices = p))
}



dfnotl <- function(dfcf,i,j){
  aux <- dfcf[-i]
  #If we are not in a triangle
  if (i != j){
    #Observations not in K=i and K=j (i.e. not in I_l)
    aux <- dfcf[-c(i,j)]
    if (length(aux) == 0){
      stop(paste("Dataframe with observations not in fold",i,
                 "and not in fold",print(j),"is empty. Consider using
                        more folds"))
    }
  }
  #ldply applies function to each element of a list and then combines results
  #in a dataframe, if no function is specified it combines everything in a
  #dataframe
  aux <- plyr::ldply(aux)
}

# weighted.mean2 <- function(X,weights = NULL){
#   if(is.null(weights)){weighted.mean(X)}
#   else{weighted.mean(X,weights)}
# }


mmdjt <- function(X,Y){
  nx <- length(X)
  ny <- length(Y)

  xx <- sapply(1:nx, function(u){
    a <- sum(exp(-0.5*(X[u] - X[-u])^2))
  })
  kxx <- sum(xx)

  yy <- sapply(1:ny, function(u){
    a <- sum(exp(-0.5*(Y[u] - Y[-u])^2))
  })
  kyy <- sum(yy)

  xy <- sapply(1:nx, function(u){
    a <- sum(exp(-0.5*(X[u] - Y)^2))
  })
  kxy <- sum(xy)

  mmd <- (1/(nx*(nx-1)))*kxx + (1/(ny*(ny-1)))*kyy - (2/(nx*ny))*kxy
}


fair_infeasible <- function(FVs,
                            A,pac1,pac0){

  n <- length(FVs)
  n1 <- n - 1
  Ppac0 <- mean(pac0)
  Ppac1 <- mean(pac1)
  fairpi <- sapply(1:n1, function(u){
    j <- u + 1
    a <- (2/(n*n1))*sum((exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u])*(pac0[j:n]))/(Ppac0*Ppac0) +
                          (exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u])*(pac1[j:n]))/(Ppac1*Ppac1) -
                          (exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac1[u])*(pac0[j:n]))/(Ppac1*Ppac0) -
                          (exp(-0.5*(FVs[u] - FVs[j:n])^2)*(pac0[u])*(pac1[j:n]))/(Ppac0*Ppac1))

  })
  fairpi <- sum(fairpi)
}
