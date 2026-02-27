#' MLU Hyperparameter Tuning
#'
#' Performs hyperparameter tuning for pairwise machine learning models across
#' multiple algorithms including Random Forest, XGBoost, CatBoost, and neural networks.
#' Uses symmetric features: levels (sum) and distances (absolute differences).
#'
#' @param X Dataframe containing all the features on which the model was estimated.
#' @param Y Vector containing the labels for which the model was estimated.
#' @param f Function of Yi and Yj defining the dependent variable.
#' @param ML Which machine learners to tune.
#' @param Kcv Number of folds for cross-validation.
#' @param rf.cf.ntree.grid How many trees should be grown when using RF or CIF.
#' @param rf.depth.grid How deep should trees be grown in RF. Inf means full depth (NULL in ranger).
#' @param mtry.grid How many variables to consider at each split in RF. Defaults floor(sqrt(ncol(X))) and floor(ncol(X)/3) are always tried.
#' @param cf.depth How deep should trees be grown in CIF (Inf is full depth, default from partykit).
#' @param ensemblefolds.grid Integer specifying how many folds to use in ensemble methods such as OLSensemble or SuperLearner.
#' @param polynomial.Lasso.grid Degree of polynomial to be fitted when using Lasso.
#' @param polynomial.Ridge.grid Degree of polynomial to be fitted when using Ridge.
#' @param polynomial.Logit_lasso.grid Degree of polynomial to be fitted when using Logit_lasso.
#' @param polynomial.OLS.grid Degree of polynomial to be fitted when using OLS.
#' @param polynomial.NLLS_exp.grid Degree of polynomial to be fitted when using NLLS_exp.
#' @param polynomial.loglin.grid Degree of polynomial to be fitted when using loglin.
#' @param xgb.nrounds.grid Number of boosting iterations for XGBoost.
#' @param xgb.eta.grid Learning rate for XGBoost.
#' @param xgb.depth.grid Maximum depth of trees for XGBoost.
#' @param xgb.child_weight.grid Minimum sum of instance weight needed in a child for XGBoost.
#' @param xgb.subsample.grid Subsample ratio of the training instances for XGBoost.
#' @param xgb.colsample.grid Subsample ratio of columns when constructing each tree for XGBoost.
#' @param cb.iterations.grid Number of boosting iterations for CatBoost.
#' @param cb.depth.grid Depth of the trees for CatBoost.
#' @param cb.learning_rate.grid Learning rate for CatBoost.
#' @param cb.l2_leaf_reg.grid L2 regularization coefficient for CatBoost.
#' @param torch.hidden_units.grid List of hidden layer configurations for neural networks.
#' @param torch.learning_rate.grid Learning rates for neural network training.
#' @param torch.dropout.grid Dropout rates for neural networks.
#' @param torch.epochs.grid Number of training epochs for neural networks.
#' @param var_penalization Penalization parameter for variance in RMSE calculation. Default is 0.
#' @param OLSensemble String vector specifying which learners should be used in OLS ensemble method.
#' @param SL.library String vector specifying which learners should be used in SuperLearner.
#' @param subsample Either NULL (use all data), a proportion between 0 and 1, or an integer number of observations to randomly subsample before creating pairs. This can significantly reduce computational time for tuning. Default is NULL.
#' @param verbose Logical specifying whether to print progress. Default is TRUE.
#'
#' @return Tuning results from ML::MLtuning().
#'
#' @export
mlutuning <- function(X, Y, f,
                       ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Torch", "Logit_lasso", "OLS", "NLLS_exp", "loglin", "OLSensemble"),
                       Kcv = 5,
                       rf.cf.ntree.grid = c(100, 300, 500),
                       rf.depth.grid = c(2, 4, 6, Inf),
                       mtry.grid = c(1, 3, 5),
                       cf.depth.grid = c(2, 4, Inf),
                       ensemblefolds.grid = c(5, 10),
                       polynomial.Lasso.grid = c(1, 2, 3),
                       polynomial.Ridge.grid = c(1, 2, 3),
                       polynomial.Logit_lasso.grid = c(1, 2, 3),
                       polynomial.OLS.grid = c(1, 2, 3),
                       polynomial.NLLS_exp.grid = c(1, 2, 3),
                       polynomial.loglin.grid = c(1, 2, 3),
                       xgb.nrounds.grid = c(100, 200, 500),
                       xgb.eta.grid = c(0.01, 0.1, 0.3),
                       xgb.depth.grid = c(1, 3, 6),
                       xgb.child_weight.grid = c(1, 5, 10),
                       xgb.subsample.grid = c(0.5, 0.75, 1),
                       xgb.colsample.grid = c(0.5, 0.75, 1),
                       cb.iterations.grid = c(100, 500, 1000),
                       cb.depth.grid = c(1, 3, 6, 10),
                       cb.learning_rate.grid = c(0.01, 0.1, 0.3),
                       cb.l2_leaf_reg.grid = c(1, 3, 5, 7, 9),
                       torch.hidden_units.grid = list(c(64, 32), c(128, 64), c(256, 128, 64)),
                       torch.learning_rate.grid = c(0.001, 0.01, 0.1),
                       torch.dropout.grid = c(0, 0.2, 0.5),
                       torch.epochs.grid = c(50, 100, 200),
                       var_penalization = 0,
                       OLSensemble = NULL,
                       SL.library = NULL,
                       subsample = NULL,
                       verbose = TRUE) {
  
  # Subsample data if requested
  n_original <- length(Y)
  if (!is.null(subsample)) {
    if (subsample > 0 && subsample < 1) {
      # Treat as proportion
      n_subsample <- floor(n_original * subsample)
    } else if (subsample >= 1 && subsample <= n_original) {
      # Treat as number of observations
      n_subsample <- floor(subsample)
    } else {
      stop("subsample must be NULL, a proportion (0-1), or an integer between 1 and n")
    }
    
    if (n_subsample < n_original) {
      set.seed(123)  # For reproducibility
      idx <- sample(1:n_original, n_subsample, replace = FALSE)
      X <- X[idx, , drop = FALSE]
      Y <- Y[idx]
      
      if (verbose) {
        cat(sprintf("\nSubsampling: %d -> %d observations (%.1f%%)\n", 
                    n_original, n_subsample, 100 * n_subsample / n_original))
      }
    }
  }
  
  # Create pairwise combinations
  n <- length(Y)
  n1 <- n - 1
  n_pairs <- n * (n - 1) * 0.5
  
  if (verbose) {
    cat("\n======================================\n")
    cat("MLU Hyperparameter Tuning\n")
    cat("======================================\n")
    if (!is.null(subsample) && n < n_original) {
      cat(sprintf("Original observations: %d (subsampled to %d)\n", n_original, n))
    } else {
      cat(sprintf("Observations: %d\n", n))
    }
    cat(sprintf("Pairs created: %d\n", n_pairs))
  }
  
  XX <- matrix(0, n_pairs, ncol(X) * 2)
  YY <- rep(0, n_pairs)
  cnt <- 0
  for (i in 1:n1) {
    j1 <- i + 1
    for (j in j1:n) {
      cnt <- cnt + 1
      Xi <- as.numeric(X[i, ])
      Xj <- as.numeric(X[j, ])
      XX[cnt, ] <- c(Xi + Xj, abs(Xi - Xj))
      YY[cnt] <- f(Y[i], Y[j])
    }
  }
  
  if (verbose) {
    cat(sprintf("Features created: %d (symmetric: levels + distances)\n", ncol(XX)))
    cat(sprintf("\nTuning ML algorithms: %s\n", paste(ML, collapse = ", ")))
    cat(sprintf("Cross-validation folds: %d\n", Kcv))
    
    # Print grid sizes
    for (ml in ML) {
      if (ml == "Lasso") {
        cat(sprintf("  - Lasso: testing %d polynomial degrees\n", length(polynomial.Lasso.grid)))
      } else if (ml == "Ridge") {
        cat(sprintf("  - Ridge: testing %d polynomial degrees\n", length(polynomial.Ridge.grid)))
      } else if (ml == "RF") {
        n_combinations <- length(rf.cf.ntree.grid) * length(rf.depth.grid) * length(mtry.grid)
        cat(sprintf("  - RF: testing %d combinations (trees × depth × mtry)\n", n_combinations))
      } else if (ml == "XGB") {
        n_combinations <- length(xgb.nrounds.grid) * length(xgb.depth.grid)
        cat(sprintf("  - XGB: testing %d combinations (rounds × depth)\n", n_combinations))
      } else if (ml == "CB") {
        n_combinations <- length(cb.iterations.grid) * length(cb.depth.grid)
        cat(sprintf("  - CB: testing %d combinations (iterations × depth)\n", n_combinations))
      } else if (ml == "Torch") {
        n_combinations <- length(torch.hidden_units.grid) * length(torch.learning_rate.grid) * 
                         length(torch.dropout.grid) * length(torch.epochs.grid)
        cat(sprintf("  - Torch: testing %d combinations (architecture × lr × dropout × epochs)\n", 
                    n_combinations))
      }
    }
    cat("======================================\n\n")
  }
  
  ML::MLtuning(as.data.frame(XX), YY, ML = ML, Kcv = Kcv,
               rf.cf.ntree.grid = rf.cf.ntree.grid,
               rf.depth.grid = rf.depth.grid,
               mtry.grid = mtry.grid,
               cf.depth.grid = cf.depth.grid,
               ensemblefolds.grid = ensemblefolds.grid,
               polynomial.Lasso.grid = polynomial.Lasso.grid,
               polynomial.Ridge.grid = polynomial.Ridge.grid,
               polynomial.Logit_lasso.grid = polynomial.Logit_lasso.grid,
               polynomial.OLS.grid = polynomial.OLS.grid,
               polynomial.NLLS_exp.grid = polynomial.NLLS_exp.grid,
               polynomial.loglin.grid = polynomial.loglin.grid,
               xgb.nrounds.grid = xgb.nrounds.grid,
               xgb.max.depth.grid = xgb.depth.grid,
               cb.iterations.grid = cb.iterations.grid,
               cb.depth.grid = cb.depth.grid,
               torch.hidden_units.grid = torch.hidden_units.grid,
               torch.lr.grid = torch.learning_rate.grid,
               torch.dropout.grid = torch.dropout.grid,
               torch.epochs.grid = torch.epochs.grid,
               var_penalization = var_penalization,
               OLSensemble = OLSensemble,
               SL.library = SL.library,
               verbose = verbose)
}
