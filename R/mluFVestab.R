#' MLU Fitted Value Estimation for All Treatment Combinations
#'
#' Estimates fitted values for pairwise ML models returning all four treatment
#' assignment combinations: (1,1), (1,0), (0,1), and (0,0) where 1 = treated
#' and 0 = untreated. Uses symmetric features: levels (sum) and distances (absolute differences).
#'
#' @param model A model object from mlumodest.
#' @param Xi Covariate dataframe or matrix for training data (observations i).
#' @param Yi Outcome vector for training data (observations i).
#' @param Xnewi Covariate dataframe or matrix for new data (observations i).
#' @param Ynewi Outcome vector for new data (observations i).
#' @param Xj Covariate dataframe or matrix for training data (observations j, for "square" mode).
#' @param Yj Outcome vector for training data (observations j, for "square" mode).
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
#' @return List with four elements:
#'   \item{fv11}{Fitted values when both units are treated}
#'   \item{fv10}{Fitted values when first unit is treated, second is untreated}
#'   \item{fv01}{Fitted values when first unit is untreated, second is treated}
#'   \item{fv00}{Fitted values when both units are untreated}
#'
#' @export
mluFVestab <- function(model,
                        Xi,
                        Yi,
                        Xnewi,
                        Ynewi,
                        Xj = NULL,
                        Yj = NULL,
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
      Xnewi <- Xnewi[idx_i, , drop = FALSE]
      Ynewi <- Ynewi[idx_i]
    }
    
    # For square mode, also subsample j
    if (shape == "square" && !is.null(Xnewj)) {
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
    XXnew11 <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    XXnew10 <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    XXnew01 <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    XXnew00 <- matrix(0, n * (n - 1) * 0.5, ncol(Xnewi) * 2 + 2)
    YYnew <- rep(0, n * (n - 1) * 0.5)
    cnt <- 0
    for (i in 1:n1) {
      j1 <- i + 1
      for (j in j1:n) {
        cnt <- cnt + 1
        Xi <- as.numeric(Xnewi[i, ])
        Xj <- as.numeric(Xnewi[j, ])
        Xsum <- Xi + Xj
        Xdiff <- abs(Xi - Xj)
        XXnew11[cnt, ] <- c(2, Xsum, 0, Xdiff)
        XXnew10[cnt, ] <- c(1, Xsum, 1, Xdiff)
        XXnew01[cnt, ] <- c(1, Xsum, 1, Xdiff)
        XXnew00[cnt, ] <- c(0, Xsum, 0, Xdiff)
        YYnew[cnt] <- f(Ynewi[i], Ynewi[j])
      }
    }
    fv11 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew11), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv10 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew10), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv01 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew01), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv00 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew00), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    return(list("fv11" = fv11, "fv10" = fv10, "fv01" = fv01, "fv00" = fv00))
  } else if (shape == "square") {
    ni <- length(Ynewi)
    nj <- length(Ynewj)
    XX <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    YY <- rep(0, ni * nj)
    XXnew11 <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    XXnew10 <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    XXnew01 <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    XXnew00 <- matrix(0, ni * nj, ncol(Xnewi) * 2 + 2)
    YYnew <- rep(0, ni * nj)
    cnt <- 0
    for (i in 1:ni) {
      for (j in 1:nj) {
        cnt <- cnt + 1
        Xi <- as.numeric(Xnewi[i, ])
        Xj <- as.numeric(Xnewj[j, ])
        Xsum <- Xi + Xj
        Xdiff <- abs(Xi - Xj)
        XXnew11[cnt, ] <- c(2, Xsum, 0, Xdiff)
        XXnew10[cnt, ] <- c(1, Xsum, 1, Xdiff)
        XXnew01[cnt, ] <- c(1, Xsum, 1, Xdiff)
        XXnew00[cnt, ] <- c(0, Xsum, 0, Xdiff)
        YYnew[cnt] <- f(Ynewi[i], Ynewj[j])
      }
    }
    fv11 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew11), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv10 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew10), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv01 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew01), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    fv00 <- ML::FVest(model, as.data.frame(XX), YY, as.data.frame(XXnew00), YYnew, ML = ML,
                      polynomial.Lasso = polynomial.Lasso,
                      polynomial.Ridge = polynomial.Ridge,
                      polynomial.Logit_lasso = polynomial.Logit_lasso,
                      polynomial.OLS = polynomial.OLS,
                      polynomial.NLLS_exp = polynomial.NLLS_exp,
                      polynomial.loglin = polynomial.loglin,
                      coefs = coefs)
    return(list("fv11" = fv11, "fv10" = fv10, "fv01" = fv01, "fv00" = fv00))
  }
}
