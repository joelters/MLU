#' MLU Model Estimation
#'
#' Estimates a machine learning model for pairwise data structures.
#' Creates all possible pairs (i,j) where i < j and fits a model using
#' symmetric features: levels (sum of features) and distances (absolute differences).
#'
#' @param D Treatment assignment vector.
#' @param X Covariate dataframe or matrix.
#' @param Y Outcome vector.
#' @param f Function of Yi and Yj defining the dependent variable for the pair.
#' @param ML String vector specifying which machine learners to use.
#'   Options: "Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Torch", "Logit_lasso", "OLS", "NLLS_exp", "loglin", "SL", "OLSensemble".
#' @param OLSensemble String vector specifying which learners should be used in OLS ensemble method.
#' @param SL.library String vector specifying which learners should be used in SuperLearner.
#' @param rf.cf.ntree How many trees should be grown when using RF or CIF. Default is 500.
#' @param rf.depth How deep should trees be grown in RF. NULL is default from ranger. Default is NULL.
#' @param mtry How many variables to consider at each split in RF. Default is max(floor(ncol(X)/3), 1).
#' @param cf.depth How deep should trees be grown in CIF. Inf is full depth. Default is Inf.
#' @param polynomial.Lasso Degree of polynomial to be fitted when using Lasso. 1 just fits the input X. 2 squares all variables and adds all pairwise interactions. Default is 1.
#' @param polynomial.Ridge Degree of polynomial to be fitted when using Ridge. See polynomial.Lasso for more info. Default is 1.
#' @param polynomial.Logit_lasso Degree of polynomial to be fitted when using Logit_lasso. See polynomial.Lasso for more info. Default is 1.
#' @param polynomial.OLS Degree of polynomial to be fitted when using OLS. See polynomial.Lasso for more info. Default is 1.
#' @param polynomial.NLLS_exp Degree of polynomial to be fitted when using NLLS_exp. See polynomial.Lasso for more info. Default is 1.
#' @param polynomial.loglin Degree of polynomial to be fitted when using loglin. See polynomial.Lasso for more info. Default is 1.
#' @param xgb.nrounds Integer specifying how many rounds to use in XGB. Default is 200.
#' @param xgb.max.depth Integer specifying how deep trees should be grown in XGB. Default is 6.
#' @param cb.iterations The maximum number of trees that can be built in CB. Default is 500.
#' @param cb.depth The depth of the trees in CB. Default is 6.
#' @param torch.epochs Integer specifying the number of epochs for Torch neural network. Default is 50.
#' @param torch.hidden_units Numeric vector specifying the number of neurons in each hidden layer of Torch neural network. Default is c(64, 32).
#' @param torch.lr Numeric value specifying the learning rate for Torch neural network. Default is 0.01.
#' @param torch.dropout Numeric value between 0 and 1 specifying the dropout rate for Torch neural network. Default is 0.2.
#' @param ensemblefolds Integer specifying how many folds to use in ensemble methods such as OLSensemble or SuperLearner. Default is 10.
#' @param start_nlls List with the starting values for NLLS_exp parameters. Default is log(mean(Y)) for intercept and zero for rest.
#' @param subsample Either NULL (use all data), a proportion between 0 and 1, or an integer number of observations to randomly subsample before creating pairs. This can significantly reduce computational time. Default is NULL.
#' @param verbose Logical specifying whether to print progress and model information. Default is TRUE.
#'
#' @details The function creates pairwise features using symmetry:
#'   For each pair (i,j), the features are constructed as:
#'   [D[i] + D[j], X[i] + X[j], |D[i] - D[j]|, |X[i] - X[j]|]
#'   This representation respects the symmetric nature of pairwise data.
#'
#' @return A model object from ML::modest().
#'
#' @export
mlumodest <- function(D, X, Y, f,
                       ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Torch", "Logit_lasso", "OLS", "NLLS_exp", "loglin", "SL", "OLSensemble"),
                       OLSensemble = NULL,
                       SL.library = NULL,
                       rf.cf.ntree = 500,
                       rf.depth = NULL,
                       mtry = NULL,
                       cf.depth = Inf,
                       polynomial.Lasso = 1,
                       polynomial.Ridge = 1,
                       polynomial.Logit_lasso = 1,
                       polynomial.OLS = 1,
                       polynomial.NLLS_exp = 1,
                       polynomial.loglin = 1,
                       xgb.nrounds = 200,
                       xgb.max.depth = 6,
                       cb.iterations = 500,
                       cb.depth = 6,
                       torch.epochs = 50,
                       torch.hidden_units = c(64, 32),
                       torch.lr = 0.01,
                       torch.dropout = 0.2,
                       ensemblefolds = 10,
                       start_nlls = NULL,
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
      D <- D[idx]
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
    cat("MLU Model Estimation\n")
    cat("======================================\n")
    if (!is.null(subsample) && n < n_original) {
      cat(sprintf("Original observations: %d (subsampled to %d)\n", n_original, n))
    } else {
      cat(sprintf("Observations: %d\n", n))
    }
    cat(sprintf("Pairs created: %d\n", n_pairs))
  }
  
  XX <- matrix(0, n_pairs, ncol(X) * 2 + 2)
  YY <- rep(0, n_pairs)
  cnt <- 0
  for (i in 1:n1) {
    j1 <- i + 1
    for (j in j1:n) {
      cnt <- cnt + 1
      Xi <- as.numeric(X[i, ])
      Xj <- as.numeric(X[j, ])
      XX[cnt, ] <- c(D[i] + D[j], Xi + Xj, abs(D[i] - D[j]), abs(Xi - Xj))
      YY[cnt] <- f(Y[i], Y[j])
    }
  }
  
  # Set default mtry if not specified
  if (is.null(mtry)) {
    mtry <- max(floor(ncol(XX)/3), 1)
  }
  
  if (verbose) {
    cat(sprintf("Features created: %d (symmetric: levels + distances)\n", ncol(XX)))
    cat(sprintf("\nFitting ML algorithm: %s\n", paste(ML, collapse = ", ")))
    
    # Print key parameters based on ML method
    if ("Lasso" %in% ML || "Ridge" %in% ML || "Logit_lasso" %in% ML || "OLS" %in% ML) {
      poly_param <- if ("Lasso" %in% ML) polynomial.Lasso 
                    else if ("Ridge" %in% ML) polynomial.Ridge
                    else if ("Logit_lasso" %in% ML) polynomial.Logit_lasso
                    else polynomial.OLS
      cat(sprintf("  - Polynomial degree: %d\n", poly_param))
    }
    if ("RF" %in% ML || "CIF" %in% ML) {
      cat(sprintf("  - Number of trees: %d\n", rf.cf.ntree))
      if (!is.null(rf.depth)) cat(sprintf("  - Tree depth: %d\n", rf.depth))
      cat(sprintf("  - mtry: %d\n", mtry))
    }
    if ("XGB" %in% ML) {
      cat(sprintf("  - Rounds: %d, Max depth: %d\n", xgb.nrounds, xgb.max.depth))
    }
    if ("CB" %in% ML) {
      cat(sprintf("  - Iterations: %d, Depth: %d\n", cb.iterations, cb.depth))
    }
    if ("Torch" %in% ML) {
      cat(sprintf("  - Hidden units: [%s]\n", paste(torch.hidden_units, collapse = ", ")))
      cat(sprintf("  - Learning rate: %.4f, Dropout: %.2f, Epochs: %d\n", 
                  torch.lr, torch.dropout, torch.epochs))
    }
    if ("SL" %in% ML) {
      cat(sprintf("  - SuperLearner library: %s\n", paste(SL.library, collapse = ", ")))
    }
    if ("OLSensemble" %in% ML) {
      cat(sprintf("  - OLS Ensemble: %s\n", paste(OLSensemble, collapse = ", ")))
      cat(sprintf("  - Ensemble folds: %d\n", ensemblefolds))
    }
    cat("======================================\n\n")
  }
  
  ML::modest(as.data.frame(XX), YY, ML = ML,
             OLSensemble = OLSensemble,
             SL.library = SL.library,
             rf.cf.ntree = rf.cf.ntree,
             rf.depth = rf.depth,
             mtry = mtry,
             cf.depth = cf.depth,
             polynomial.Lasso = polynomial.Lasso,
             polynomial.Ridge = polynomial.Ridge,
             polynomial.Logit_lasso = polynomial.Logit_lasso,
             polynomial.OLS = polynomial.OLS,
             polynomial.NLLS_exp = polynomial.NLLS_exp,
             polynomial.loglin = polynomial.loglin,
             xgb.nrounds = xgb.nrounds,
             xgb.max.depth = xgb.max.depth,
             cb.iterations = cb.iterations,
             cb.depth = cb.depth,
             torch.epochs = torch.epochs,
             torch.hidden_units = torch.hidden_units,
             torch.lr = torch.lr,
             torch.dropout = torch.dropout,
             ensemblefolds = ensemblefolds,
             start_nlls = start_nlls)
}
