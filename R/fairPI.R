fairPI <- function(Y,
                   X,
                   A,
                   ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "Logit_lasso"),
                   fairness = c("SP", "EO", "PPE"),
                   relative = TRUE,
                   rf.ntree = 500,
                   polynomial = 1,
                   rf.depth = NULL){

  # select additional "nuisance" according to fairness definition
  if (fairness == "SP"){
    pac0 <- A == 0
    pac1 <- A == 1
  }
  else if (fairness == "EO"){
    pac0 <- (A == 0)*(Y == 1)
    pac1 <- (A == 1)*(Y == 1)
  }
  else if (fairness == "PPE"){
    pac0 <- (A == 0)*(Y == 0)
    pac1 <- (A == 1)*(Y == 0)
  }

  fair <- mlfair(X,
                 Y,
                 A,
                 pac1 = pac1,
                 pac0 = pac0,
                 relative = relative,
                 ML = ML,
                 fitted_values = TRUE,
                 rf.ntree = rf.ntree,
                 rf.depth = rf.depth)
  fairpi <- fair$fairMMD

  FVs <- fair$FVs


  if (relative == TRUE){
    Dmax <- 2 - 2*exp(-0.5)
    fair_res <- as.matrix(c(fairpi[1],fairpi[2]))
    rownames(fair_res) <- c("Fair", "Fair_rel")
  }
  else{
    fair_res <- as.matrix(c(fairpi))
    rownames(fair_res) <- c("Fair")
  }
  return(data.frame(FairMMD = fair_res))
}
