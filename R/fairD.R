fairD <- function(Y,
                  X,
                  A,
                  CFit = TRUE,
                  npart = 5,
                  ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "Logit_lasso"),
                  ensemble = c("SL.Lasso","SL.Ridge","SL.RF","SL.CIF",
                               "SL.XGB","SL.Logit_lasso","SL.CB"),
                  fairness = c("SP", "EO", "PPE"),
                  relative = TRUE,
                  sterr = TRUE,
                  fitted_values = TRUE,
                  parallel = FALSE,
                  rf.ntree = 500,
                  rf.depth = NULL,
                  polynomial = 1,
                  verbose = FALSE){

  #method takes the value of RF_type
  ML = match.arg(ML)
  na0 <- sum(A == 0)
  na1 <- sum(A == 1)
  if (fairness == "SP"){
    Ppac0 <- mean(A == 0)
    Ppac1 <- mean(A == 1)
  }
  else if (fairness == "EO"){
    Ppac0 <- mean((A == 0)*(Y == 1))
    Ppac1 <- mean((A == 1)*(Y == 1))
  }
  else if (fairness == "PPE"){
    Ppac0 <- mean((A == 0)*(Y == 0))
    Ppac1 <- mean((A == 1)*(Y == 0))
  }
  n <- length(Y)
  if (CFit == FALSE){
    #Model and FVs estimation
    m <- ML::MLest(X,
                   Y,
                   ML,
                   ensemble = ensemble,
                   rf.cf.ntree = rf.ntree,
                   rf.depth = rf.depth,
                   polynomial = polynomial,
                   FVs = TRUE)
    model <- m$model
    FVs <- m$FVs
    if(sum(FVs < 0 | FVs > 1) != 0){
      warning("There are estimated probabilities outside unit interval")
    }

    # select additional "nuisance" according to fairness definition
    if (fairness == "SP"){
      pac0 <- A == 0
      pac1 <- A == 1
      pacx1 <- pac1
      pacx0 <- pac0
    }
    else if (fairness == "EO"){
      pac0 <- (A == 0)*(Y == 1)
      pac1 <- (A == 1)*(Y == 1)
      pacx0 <- (A == 0)*FVs
      pacx1 <- (A == 1)*FVs
    }
    else if (fairness == "PPE"){
      pac0 <- (A == 0)*(Y == 0)
      pac1 <- (A == 1)*(Y == 0)
      pacx0 <- (A == 0)*(1 - FVs)
      pacx1 <- (A == 1)*(1 - FVs)
    }

    #Debiased Fairness
    # browser()
    fairdeb <- fair_deb(Y, FVs, pac1, pac0, pacx1, pacx0, A)
    th11 <- fairdeb['th11']
    th00 <- fairdeb['th00']
    th10 <- fairdeb['th10']
    th01 <- fairdeb['th01']
    fairdeb <- th00 + th11 - th10 - th01
    #SEs
    if (sterr == TRUE){
      se <- se_deb(Y, FVs, A, pac1, pac0, pacx1, pacx0,
                    th11, th00, th10, th01)
    } else{se <- NULL}
    return(data.frame(fairMMD = fairdeb, se = se, row.names = ""))
  }
  else{
    nn <- length(Y)
    df <- dplyr::as_tibble(cbind(cbind(Y = Y, A = A),X))
    dfcf <- SP(df,npart)
    dfcf <- dfcf$dfsp
    numcf <- 0
    # SM <- vector(mode = "list", npart)
    count <- 0
    cnt <- 0
    RMSE1 <- rep(0,npart)
    CFind <- NULL
    for (i in 1:npart){
      count <- count + 1
      j1 <- i
      for (j in j1:npart){
        cnt <- cnt + 1
        CFind[[cnt]] <- c(i,j,cnt,count)
      }
    }
    if (parallel == FALSE){
      rescf <- lapply(CFind,function(u){
        i <- u[1]
        j <- u[2]
        cnt <- u[3]
        count <- u[4]

        #Create dataframe with all observations not in I_l
        aux <- dfnotl(dfcf,i,j)

        #Train model with observations not in I_l (not in K = i or K = j)
        m <- ML::MLest(dplyr::select(aux,-c(Y,A)),
                       aux$Y,
                       ML,
                       ensemble = ensemble,
                       rf.cf.ntree = rf.ntree,
                       rf.depth = rf.depth,
                       polynomial = polynomial,
                       FVs = FALSE)
        model <- m$model

        #Predict fitted values for obs in Ci and Cj
        #using model trained with obs not in Ci or Cj
        #Y1 are income observations in Ci, FVs1 are
        #predictions on Ci with model trained with obs
        #not in Ci or Cj
        #Y2 are income observations in Cj, FVs2 are
        #predictions on Cj with model trained with obs
        #not in Ci or Cj

        #i <= j (we want to sum across all i,j, but we are
        #doing i <= j in the loop, we are going to do
        #i > j manually whenever in the loop i!=j by just
        #switching the roles of i and j, this way we dont
        #train the model twice (once for (i,j) and another for (j,i)))
        X1 <- dplyr::select(dfcf[[i]], -c(Y,A))
        Y1 <- dfcf[[i]]$Y
        A1 <- dfcf[[i]]$A
        FVs1 <- ML::FVest(model,X,Y,X1,Y1,ML)

        # select additional "nuisance" according to fairness definition
        if (fairness == "SP"){
          pac1_0 <- A1 == 0
          pac1_1 <- A1 == 1
          pacx1_0 <- pac1_0
          pacx1_1 <- pac1_1
        }
        else if (fairness == "EO"){
          pac1_0 <- (A1 == 0)*(Y1 == 1)
          pac1_1 <- (A1 == 1)*(Y1 == 1)
          pacx1_0 <- (A1 == 0)*FVs1
          pacx1_1 <- (A1 == 1)*FVs1
        }
        else if (fairness == "PPE"){
          pac1_0 <- (A1 == 0)*(Y1 == 0)
          pac1_1 <- (A1 == 1)*(Y1 == 0)
          pacx1_0 <- (A1 == 0)*(1 - FVs1)
          pacx1_1 <- (A1 == 1)*(1 - FVs1)
        }

        #If we are in a square
        if (i != j){
          Y2 <- dfcf[[j]]$Y
          A2 <- dfcf[[j]]$A
          X2 <- dplyr::select(dfcf[[j]], -c(Y,A))
          FVs2 <- ML::FVest(model,X,Y,X2,Y2,ML)
          if(sum(FVs1 < 0 | FVs1 > 1) != 0 |
             sum(FVs2 < 0 | FVs2 > 1) != 0){
            warning("There are estimated probabilities outside unit interval")
          }

          # select additional "nuisance" according to fairness definition
          if (fairness == "SP"){
            pac2_0 <- A2 == 0
            pac2_1 <- A2 == 1
            pacx2_0 <- pac2_0
            pacx2_1 <- pac2_1
          }
          else if (fairness == "EO"){
            pac2_0 <- (A2 == 0)*(Y2 == 1)
            pac2_1 <- (A2 == 1)*(Y2 == 1)
            pacx2_0 <- (A2 == 0)*FVs2
            pacx2_1 <- (A2 == 1)*FVs2
          }
          else if (fairness == "PPE"){
            pac2_0 <- (A2 == 0)*(Y2 == 0)
            pac2_1 <- (A2 == 1)*(Y2 == 0)
            pacx2_0 <- (A2 == 0)*(1 - FVs2)
            pacx2_1 <- (A2 == 1)*(1 - FVs2)
          }

          b <- fairdnumsq(Y1, Y2, FVs1, FVs2, A1, A2,
                          pac1_1, pac1_0, pacx1_1, pacx1_0,
                          pac2_1, pac2_0, pacx2_1, pacx2_0,
                          Ppac0, Ppac1)
        }

        #In Triangle
        fvss <- rep(0,length(Y1))
        RMSE1 <- 0
        if (i == j){
          # b <- fairdnumtr2(Y1, FVs1, A1, na0, na1)
          # b <- fairdnumtr3(Y1, FVs1, A1,
          #                  pac1_1, pac1_0, pacx1_1, pacx1_0,
          #                  na0, na1)
          b <- fairdnumtr(Y1, FVs1, A1,
                           pac1_1, pac1_0, pacx1_1, pacx1_0,
                           Ppac0, Ppac1)
          RMSE1 <- (length(Y1)/length(Y))*
            sqrt(mean(((Y1 - FVs1)^2)))
          fvss <- FVs1
          if (verbose == TRUE){
            print(paste(round(100*cnt/(npart*(npart + 1)/2),2),"% completed"))
          }
        }
        #Add the term in the numerator of the estimator
        #corresponding to this block
        # return(list("fairrmse" = data.frame("fair" = b,
        #                                     "RMSE1" = RMSE1),
        #             FVs = fvss))
        return(list("fairrmse" = data.frame("th00" = b["th00"],
                                            "th11" = b["th11"],
                                            "th10" = b["th10"],
                                            "th01" = b["th01"],
                                            "RMSE1" = RMSE1),
                    FVs = fvss))
      })
    }
    else if (parallel == TRUE){
      n.cores <- parallel::detectCores()
      clust <- parallel::makeCluster(n.cores)
      parallel::clusterEvalQ(clust, set.seed(123))
      parallel::clusterExport(clust, c("dfcf","npart","rf.cf.ntree",
                                       "rf.depth","polynomial", "Ppac1",
                                       "Ppac0"),
                              envir=environment())
      rescf <- parallel::parLapply(clust, CFind, function(u){
        i <- u[1]
        j <- u[2]
        cnt <- u[3]
        count <- u[4]

        #Create dataframe with all observations not in I_l
        aux <- dfnotl(dfcf,i,j)

        #Train model with observations not in I_l (not in K = i or K = j)
        m <- ML::MLest(dplyr::select(aux,-c(Y,A)),
                       aux$Y,
                       ML,
                       ensemble = ensemble,
                       rf.cf.ntree = rf.ntree,
                       rf.depth = rf.depth,
                       polynomial = polynomial,
                       FVs = FALSE)
        model <- m$model

        #Predict fitted values for obs in Ci and Cj
        #using model trained with obs not in Ci or Cj
        #Y1 are income observations in Ci, FVs1 are
        #predictions on Ci with model trained with obs
        #not in Ci or Cj
        #Y2 are income observations in Cj, FVs2 are
        #predictions on Cj with model trained with obs
        #not in Ci or Cj

        #i <= j (we want to sum across all i,j, but we are
        #doing i <= j in the loop, we are going to do
        #i > j manually whenever in the loop i!=j by just
        #switching the roles of i and j, this way we dont
        #train the model twice (once for (i,j) and another for (j,i)))
        X1 <- dplyr::select(dfcf[[i]], -c(Y,A))
        Y1 <- dfcf[[i]]$Y
        A1 <- dfcf[[i]]$A
        FVs1 <- ML::FVest(model,X,Y,X1,Y1,ML)

        # select additional "nuisance" according to fairness definition
        if (fairness == "SP"){
          pac1_0 <- A1 == 0
          pac1_1 <- A1 == 1
          pacx1_0 <- pac1_0
          pacx1_1 <- pac1_1
        }
        else if (fairness == "EO"){
          pac1_0 <- (A1 == 0)*(Y1 == 1)
          pac1_1 <- (A1 == 1)*(Y1 == 1)
          pacx1_0 <- (A1 == 0)*FVs1
          pacx1_1 <- (A1 == 1)*FVs1
        }
        else if (fairness == "PPE"){
          pac1_0 <- (A1 == 0)*(Y1 == 0)
          pac1_1 <- (A1 == 1)*(Y1 == 0)
          pacx1_0 <- (A1 == 0)*(1 - FVs1)
          pacx1_1 <- (A1 == 1)*(1 - FVs1)
        }

        #If we are in a square
        if (i != j){
          Y2 <- dfcf[[j]]$Y
          A2 <- dfcf[[j]]$A
          X2 <- dplyr::select(dfcf[[j]], -c(Y,A))
          FVs2 <- ML::FVest(model,X,Y,X2,Y2,ML)


          if(sum(FVs1 < 0 | FVs1 > 1) != 0){
            warning(paste(sum(FVs <= 0),"There are estimated probabilities outside unit interval"))
          }
          if(sum(FVs2 < 0 | FVs2 > 1) != 0){
            warning(paste(sum(FVs <= 0),"There are estimated probabilities outside unit interval"))
          }
          # select additional "nuisance" according to fairness definition
          if (fairness == "SP"){
            pac2_0 <- A2 == 0
            pac2_1 <- A2 == 1
            pacx2_0 <- pac2_0
            pacx2_1 <- pac2_1
          }
          else if (fairness == "EO"){
            pac2_0 <- (A2 == 0)*(Y2 == 1)
            pac2_1 <- (A2 == 1)*(Y2 == 1)
            pacx2_0 <- (A2 == 0)*FVs2
            pacx2_1 <- (A2 == 1)*FVs2
          }
          else if (fairness == "PPE"){
            pac2_0 <- (A2 == 0)*(Y2 == 0)
            pac2_1 <- (A2 == 1)*(Y2 == 0)
            pacx2_0 <- (A2 == 0)*(1 - FVs2)
            pacx2_1 <- (A2 == 1)*(1 - FVs2)
          }

          # b <- fairdnumsq2(Y1, Y2, FVs1, FVs2, A1, A2, na0, na1)

          b <- fairdnumsq(Y1, Y2, FVs1, FVs2, A1, A2,
                          pac1_1, pac1_0, pacx1_1, pacx1_0,
                          pac2_1, pac2_0, pacx2_1, pacx2_0,
                          Ppac0, Ppac1)
        }

        #In Triangle
        fvss <- rep(0,length(Y1))
        RMSE1 <- 0
        if (i == j){
          # b <- fairdnumtr2(Y1, FVs1, A1, na0, na1)
          b <- fairdnumtr(Y1, FVs1, A1,
                          pac1_1, pac1_0, pacx1_1, pacx1_0,
                          Ppac0, Ppac1)
          RMSE1 <- (length(Y1)/length(Y))*
            sqrt(mean((Y1 - FVs1)^2))
          fvss <- FVs1
          if (verbose == TRUE){
            print(paste(round(100*cnt/(npart*(npart + 1)/2),2),"% completed"))
          }
        }

        #Add the term in the numerator of the estimator
        #corresponding to this block
        return(list("fairrmse" = data.frame("fair" = b,
                                            "RMSE1" = RMSE1),
                    FVs = fvss))

      })
      parallel::stopCluster(clust)
    }

    rescf2 <- lapply(1:length(rescf), function(u) rescf[[u]][[1]])
    rescf3 <- do.call(rbind,rescf2)
    iii <- which(rescf3$RMSE1 != 0)
    rescf4 <- lapply(1:length(rescf), function(u) rescf[[u]][[2]])
    FVs <- unlist(rescf4[iii])

    # fairmm <- 2*sum(rescf3$fair)
    th11 <- (2/(n*(n-1)))*sum(rescf3$th11)
    th00 <- (2/(n*(n-1)))*sum(rescf3$th00)
    th10 <- (2/(n*(n-1)))*sum(rescf3$th10)
    th01 <- (2/(n*(n-1)))*sum(rescf3$th01)
    fairmm <- th00 + th11 - th10 - th01
    RMSE1 <- sum(rescf3$RMSE1)

    #FVs
    if (fitted_values == TRUE | sterr == TRUE){
      m <- ML::MLest(X, Y, ML, ensemble = ensemble, FVs = TRUE,
                     rf.cf.ntree = rf.ntree,
                     rf.depth = rf.depth,
                     polynomial = polynomial)
      FVres <- m$FVs
    } else if (fitted_values == FALSE & sterr == FALSE){
      FVres <- NULL
    }
    #SE
    if (sterr == TRUE){
      if (verbose == TRUE){
        print("Computing standard error")
      }
      # select additional "nuisance" according to fairness definition
      if (fairness == "SP"){
        pac0 <- A == 0
        pac1 <- A == 1
        pacx1 <- pac1
        pacx0 <- pac0
      }
      else if (fairness == "EO"){
        pac0 <- (A == 0)*(Y == 1)
        pac1 <- (A == 1)*(Y == 1)
        pacx0 <- (A == 0)*FVs
        pacx1 <- (A == 1)*FVs
      }
      else if (fairness == "PPE"){
        pac0 <- (A == 0)*(Y == 0)
        pac1 <- (A == 1)*(Y == 0)
        pacx0 <- (A == 0)*(1 - FVs)
        pacx1 <- (A == 1)*(1 - FVs)
      }
      # se <- se_deb(Y, FVs, A, pac1, pac0, pacx1, pacx0, fairmm)
      se <- se_deb(Y, FVs, A, pac1, pac0, pacx1, pacx0,
                    th11, th00, th10, th01)
    } else{
      se <- NULL
    }

    #Relative
    #IOp relative
    if (relative == TRUE){
      Dmax <- 2 - 2*exp(-0.5)
      fair_rel <- fairmm/Dmax
      if (sterr == TRUE){
        se_rel_fair <- se/Dmax
      } else{se_rel_fair <- NULL}
    } else{
      fair_rel <- NULL
      se_rel_fair <- NULL
    }


    jt <- rbind("Fair" = fairmm, "se" = se)
    jt_rel <- rbind("Fair_rel" = fair_rel, "se" = se_rel_fair)
    colnames(jt) <- "Fair"
    if (fitted_values == TRUE){
      return(list(fairMMD = jt, RMSE1 = RMSE1, fair_rel = jt_rel, FVs = FVres))
    } else{
      return(list(fairMMD = jt, RMSE1 = RMSE1, fair_rel = jt_rel))
    }
  }
}
