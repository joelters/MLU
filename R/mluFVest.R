#' MLU Fitted Value Estimation
#'
#' Estimates fitted values for pairwise ML models with treatment assignment.
#' Supports two modes: "triangle" for within-group pairs and "square" for
#' between-group pairs. Uses symmetric features: levels (sum) and distances (absolute differences).
#'
#' @param model A model object from mlumodest.
#' @param Di Treatment assignment vector for training data (observations i).
#' @param Xi Covariate dataframe or matrix for training data (observations i).
#' @param Yi Outcome vector for training data (observations i).
#' @param Dnewi Treatment assignment vector for new data (observations i).
#' @param Xnewi Covariate dataframe or matrix for new data (observations i).
#' @param Ynewi Outcome vector for new data (observations i).
#' @param Dj Treatment assignment vector for training data (observations j, for "square" mode).
#' @param Xj Covariate dataframe or matrix for training data (observations j, for "square" mode).
#' @param Yj Outcome vector for training data (observations j, for "square" mode).
#' @param Dnewj Treatment assignment vector for new data (observations j, for "square" mode).
#' @param Xnewj Covariate dataframe or matrix for new data (observations j, for "square" mode).
#' @param Ynewj Outcome vector for new data (observations j, for "square" mode).
#' @param f Function of Yi and Yj defining the dependent variable.
#' @param ML String vector specifying which machine learners to use.
#' @param shape Either "triangle" (within-sample pairs: i < j from same dataset) or
#'   "square" (between-sample pairs: all i,j combinations from two datasets).
#' @param polynomial.Lasso Degree of polynomial to be fitted when using Lasso. Default is 1.
#' @param polynomial.Ridge Degree of polynomial to be fitted when using Ridge. Default is 1.
#' @param polynomial.Logit_lasso Degree of polynomial to be fitted when using Logit_lasso. Default is 1.
#' @param polynomial.OLS Degree of polynomial to be fitted when using OLS. Default is 1.
#' @param polynomial.NLLS_exp Degree of polynomial to be fitted when using NLLS_exp. Default is 1.
#' @param polynomial.loglin Degree of polynomial to be fitted when using loglin. Default is 1.
#' @param coefs Optimal coefficients for OLSensemble, computed in mlumodest. Default is NULL.
#' @param subsample Either NULL (use all data), a proportion between 0 and 1, or an integer number of observations to randomly subsample before creating pairs. Applies to both i and j observations separately. Default is NULL.
#'
#' @return Fitted values from ML::FVest().
#'
#' @export
mluFVest <- function(model,
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
                      ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Torch", "Logit_lasso", "OLS", "NLLS_exp", "loglin", "SL", "OLSensemble"),
                      shape = c("triangle", "square"),
                      polynomial.Lasso = 1,
                      polynomial.Ridge = 1,
                      polynomial.Logit_lasso = 1,
                      polynomial.OLS = 1,
                      polynomial.NLLS_exp = 1,
                      polynomial.loglin = 1,
                      coefs = NULL,
                      subsample = NULL) {
  
  # Subsample data if requested
  if (!is.null(subsample)) {
    n_original_i <- length(Ynewi)
    if (subsample > 0 && subsample < 1) {
      n_subsample_i <- floor(n_original_i * subsample)
    } else if (subsample >= 1 && subsample <= n_original_i) {
      n_subsample_i <- floor(subsample)
    } else {
      stop("subsample must be NULL, a proportion (0-1), or an integer between 1 and n")
    }
    
    if (n_subsample_i < n_original_i) {
      set.seed(123)
      idx_i <- sample(1:n_original_i, n_subsample_i, replace = FALSE)
      Dnewi <- Dnewi[idx_i]
      Xnewi <- Xnewi[idx_i, , drop = FALSE]
      Ynewi <- Ynewi[idx_i]
    }
    
    # For square mode, also subsample j
    if (shape == "square" && !is.null(Dnewj)) {
      n_original_j <- length(Ynewj)
      if (subsample > 0 && subsample < 1) {
        n_subsample_j <- floor(n_original_j * subsample)
      } else if (subsample >= 1 && subsample <= n_original_j) {
        n_subsample_j <- floor(subsample)
      } else {
        n_subsample_j <- n_original_j
      }
      
      if (n_subsample_j < n_original_j) {
        set.seed(456)
        idx_j <- sample(1:n_original_j, n_subsample_j, replace = FALSE)
        Dnewj <- Dnewj[idx_j]
        Xnewj <- Xnewj[idx_j, , drop = FALSE]
        Ynewj <- Ynewj[idx_j]
      }
    }
  }
  
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
        Xi <- as.numeric(Xnewi[i, ])
        Xj <- as.numeric(Xnewi[j, ])
        XXnew[cnt, ] <- c(Dnewi[i] + Dnewi[j], Xi + Xj, abs(Dnewi[i] - Dnewi[j]), abs(Xi - Xj))
        YYnew[cnt] <- f(Ynewi[i], Ynewi[j])
      }
    }
    ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew), YYnew, ML = ML,
              polynomial.Lasso = polynomial.Lasso,
              polynomial.Ridge = polynomial.Ridge,
              polynomial.Logit_lasso = polynomial.Logit_lasso,
              polynomial.OLS = polynomial.OLS,
              polynomial.NLLS_exp = polynomial.NLLS_exp,
              polynomial.loglin = polynomial.loglin,
              coefs = coefs)
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
        Xi <- as.numeric(Xnewi[i, ])
        Xj <- as.numeric(Xnewj[j, ])
        XXnew[cnt, ] <- c(Dnewi[i] + Dnewj[j], Xi + Xj, abs(Dnewi[i] - Dnewj[j]), abs(Xi - Xj))
        YYnew[cnt] <- f(Ynewi[i], Ynewj[j])
      }
    }
    ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew), YYnew, ML = ML,
              polynomial.Lasso = polynomial.Lasso,
              polynomial.Ridge = polynomial.Ridge,
              polynomial.Logit_lasso = polynomial.Logit_lasso,
              polynomial.OLS = polynomial.OLS,
              polynomial.NLLS_exp = polynomial.NLLS_exp,
              polynomial.loglin = polynomial.loglin,
              coefs = coefs)
  }
}
