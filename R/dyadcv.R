#' Dyadic ML Cross-Validation
#'
#' Performs cross-validation for dyadic machine learning models to select
#' the best performing learner based on RMSE.
#'
#' @param X Dataframe containing all the features on which the model was estimated.
#' @param Y Vector containing the labels for which the model was estimated.
#' @param f Function of Yi and Yj defining the dependent variable.
#' @param ML String vector specifying which machine learners to use.
#' @param Kcv Number of folds for cross-validation. Default is 5.
#'
#' @return List containing ML attaining minimum RMSE and RMSE value.
#'
#' @export
dyadcv <- function(X, Y, f,
                   ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Logit_lasso"),
                   Kcv = 5) {
  # Create dyads
  n <- length(Y)
  n1 <- n - 1
  XX <- matrix(0, n * (n - 1) * 0.5, ncol(X) * 2)
  YY <- rep(0, n * (n - 1) * 0.5)
  cnt <- 0
  for (i in 1:n1) {
    j1 <- i + 1
    for (j in j1:n) {
      cnt <- cnt + 1
      XX[cnt, ] <- c(as.numeric(X[i, ]), as.numeric(X[j, ]))
      YY[cnt] <- f(Y[i], Y[j])
    }
  }
  ML::MLcv(as.data.frame(XX), YY, ML = ML, Kcv = Kcv)
}
