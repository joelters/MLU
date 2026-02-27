#' Dyadic Fitted Value Estimation
#'
#' Estimates fitted values for dyadic models with treatment assignment.
#' Supports two modes: "triangle" for within-group pairs and "square" for
#' between-group pairs.
#'
#' @param model A model object from dyadmodest.
#' @param Di Treatment assignment vector for training data (group i).
#' @param Xi Covariate dataframe or matrix for training data (group i).
#' @param Yi Outcome vector for training data (group i).
#' @param Dnewi Treatment assignment vector for new data (group i).
#' @param Xnewi Covariate dataframe or matrix for new data (group i).
#' @param Ynewi Outcome vector for new data (group i).
#' @param Dj Treatment assignment vector for training data (group j, for "square" mode).
#' @param Xj Covariate dataframe or matrix for training data (group j, for "square" mode).
#' @param Yj Outcome vector for training data (group j, for "square" mode).
#' @param Dnewj Treatment assignment vector for new data (group j, for "square" mode).
#' @param Xnewj Covariate dataframe or matrix for new data (group j, for "square" mode).
#' @param Ynewj Outcome vector for new data (group j, for "square" mode).
#' @param f Function of Yi and Yj defining the dependent variable.
#' @param ML String vector specifying which machine learners to use.
#' @param shape Either "triangle" (within-group pairs: i < j from same group) or
#'   "square" (between-group pairs: all i,j combinations from two groups).
#'
#' @return Fitted values from ML::FVest().
#'
#' @export
dyadFVest <- function(model,
                      Di,
                      Xi,
                      Yi,
                      Dnewi,
                      Xnewi,
                      Ynewi,
                      Dj = NULL,
                      Xj = NULL,
                      Yj = NULL,
                      Dnewj = NULL,
                      Xnewj = NULL,
                      Ynewj = NULL,
                      f,
                      ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Logit_lasso", "SL"),
                      shape = c("triangle", "square")) {
  if (shape == "triangle") {
    n <- length(Ynewi)
    n1 <- n - 1
    XX <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    YY <- rep(0, n * (n - 1) * 0.5)
    XXnew <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    YYnew <- rep(0, n * (n - 1) * 0.5)
    cnt <- 0
    for (i in 1:n1) {
      j1 <- i + 1
      for (j in j1:n) {
        cnt <- cnt + 1
        XXnew[cnt, ] <- c(Dnewi[i], as.numeric(Xnewi[i, ]), Dnewi[j], as.numeric(Xnewi[j, ]))
        YYnew[cnt] <- f(Ynewi[i], Ynewi[j])
      }
    }
    ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew), YYnew, ML = ML)
  } else if (shape == "square") {
    ni <- length(Ynewi)
    nj <- length(Ynewj)
    XX <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    YY <- rep(0, ni * nj)
    XXnew <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    YYnew <- rep(0, ni * nj)
    cnt <- 0
    for (i in 1:ni) {
      for (j in 1:nj) {
        cnt <- cnt + 1
        XXnew[cnt, ] <- c(Dnewi[i], as.numeric(Xnewi[i, ]), Dnewj[j], as.numeric(Xnewj[j, ]))
        YYnew[cnt] <- f(Ynewi[i], Ynewj[j])
      }
    }
    ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew), YYnew, ML = ML)
  }
}
