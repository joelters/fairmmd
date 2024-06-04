#' Estimate Maximum Mean Discrepancy (MMD)
#'
#' `fairmmd` estimates MMD for fairness in different ways. Plug in estimators are available
#' where the fitted values are estimated
#' by any machine learner among Lasso, Ridge,
#' Random Forest, Conditional Inference Forest, Extreme Gradient
#' Boosting, Catboosting. Debiased estimates based on the aforementioned
#' machine learners can also be computed with this function. In this case the
#' (optional) standard errors are derived analytically and backed by inferential
#' theory.
#'
#' @param Y is a vector containing the outcome of interest
#' @param X is a dataframe containing observable characteristics
#' @param A is a vector containing the sensitive feature
#' @param est_method whether to compute plug-in or debiased estimator
#' @param ML is a string specifying which machine learner to use
#' @param fairness vector containing three different fairness definitions: SP (
#' Statistical parity), EO (Equality of Opportunity) and PPE (Probabilistic Predictive
#' Equality).
#' @param sterr logical indicating whether standard errors should be computed
#' @param CFit logical indicating whether Cross-Fitting should be done in
#' the debiased estimators (no inferential guarantee can be given yet if FALSE)
#' @param npart in how many parts should the data be split for cross-fitting
#' @param parallel whether to parallelize cross-fitting procedure
#' @param relative logical indicating whether relative (to maximum possible MMD)
#' MMD should be computed
#' @param rf.ntree number of trees to grow if using Random Forests or Conditional
#' Inference Forests
#' @param rf.depth how deep should the Random Forest trees be. See ranger documentation.
#' @param polynomial degree of polynomial to be fitted when using Lasso, Ridge
#' or Logit Lasso. 1 just fits the input X. 2 squares all variables and adds
#' all pairwise interactions. 3 squares and cubes all variables and adds all
#' pairwise and threewise interactions...
#' @param verbose whether some progress should be reported while the function runs.
#' @returns list containing MMD estimates, RMSE of the first stage (for Debiased
#' estimates), relative MMD (if desired) and fitted values (if desired).
#' @examples
#' compas <- dplyr::select(compas, -c("probability","predicted"))
#' compas$Two_yr_Recidivism <- as.numeric(compas$Two_yr_Recidivism) - 1
#' compas$black <- compas$ethnicity == "African_American"
#'
#' mmdpi <- fairmmd(Y = compas$Two_yr_Recidivism,
#'                  X = dplyr::select(compas,-c("Two_yr_Recidivism","ethnicity")),
#'                  A = compas$black,
#'                  est_method = "Plugin",
#'                  ML = "RF",
#'                  fairness = "SP",
#'                  relative = TRUE,
#'                  sterr = FALSE)
#'
#' mmddeb <- fairmmd(Y = compas$Two_yr_Recidivism,
#'                   X = dplyr::select(compas,-c("Two_yr_Recidivism","ethnicity")),
#'                   A = compas$black,
#'                   est_method = "Debiased",
#'                   ML = "RF",
#'                   fairness = "SP",
#'                   relative = TRUE,
#'                   sterr = FALSE)
#'
#'
#' @references Escanciano, J. C., & Terschuur, J. R. (2022).
#' Debiased Semiparametric U-Statistics: Machine Learning Inference
#' on Inequality of Opportunity. arXiv preprint arXiv:2206.05235.
#' @export
fairmmd <- function(Y,
                    X,
                    A,
                    est_method = c("Plugin","Debiased"),
                    ML = c("Lasso", "Ridge", "RF", "CIF", "XGB","Logit_lasso"),
                    fairness = c("SP", "EO", "PPE"),
                    relative = TRUE,
                    sterr = TRUE,
                    CFit = TRUE,
                    npart = 5,
                    parallel = FALSE,
                    rf.ntree = 500,
                    rf.depth = NULL,
                    polynomial = 1,
                    verbose = FALSE){
  fairness = match.arg(fairness)
  if (est_method == "Plugin"){
    fair <- fairPI(Y,
                   X,
                   A,
                   ML = ML,
                   fairness = fairness,
                   relative = relative,
                   rf.ntree = rf.ntree,
                   rf.depth = rf.depth,
                   polynomial = polynomial)
  }
  else if (est_method == "Debiased"){
    fair <- fairD(Y,
                 X,
                 A,
                 CFit = CFit,
                 npart = npart,
                 ML = ML,
                 fairness = fairness,
                 relative = relative,
                 sterr = sterr,
                 parallel = parallel,
                 rf.ntree = rf.ntree,
                 rf.depth = rf.depth,
                 polynomial = polynomial,
                 verbose = verbose)
  }
}
