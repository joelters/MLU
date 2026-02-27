#' Dyadic ML Hyperparameter Tuning
#'
#' Performs hyperparameter tuning for dyadic machine learning models across
#' multiple algorithms including Random Forest, XGBoost, CatBoost, and neural networks.
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
#'
#' @return Tuning results from ML::MLtuning().
#'
#' @export
dyadtuning <- function(X, Y, f,
                       ML = c("Lasso", "Ridge", "RF", "CIF", "XGB", "CB", "Logit_lasso"),
                       Kcv = 5,
                       rf.cf.ntree.grid = c(100, 300, 500),
                       rf.depth.grid = c(2, 4, 6, Inf),
                       mtry.grid = c(1, 3, 5),
                       cf.depth.grid = c(2, 4, Inf),
                       ensemblefolds.grid = c(5, 10),
                       polynomial.Lasso.grid = c(1, 2, 3),
                       polynomial.Ridge.grid = c(1, 2, 3),
                       polynomial.Logit_lasso.grid = c(1, 2, 3),
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
                       torch.epochs.grid = c(50, 100, 200)) {
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
  
  ML::MLtuning(as.data.frame(XX), YY, ML = ML, Kcv = Kcv,
               rf.cf.ntree.grid = rf.cf.ntree.grid,
               rf.depth.grid = rf.depth.grid,
               mtry.grid = mtry.grid,
               cf.depth.grid = cf.depth.grid,
               ensemblefolds.grid = ensemblefolds.grid,
               polynomial.Lasso.grid = polynomial.Lasso.grid,
               polynomial.Ridge.grid = polynomial.Ridge.grid,
               polynomial.Logit_lasso.grid = polynomial.Logit_lasso.grid,
               xgb.nrounds.grid = xgb.nrounds.grid,
               xgb.eta.grid = xgb.eta.grid,
               xgb.depth.grid = xgb.depth.grid,
               xgb.child_weight.grid = xgb.child_weight.grid,
               xgb.subsample.grid = xgb.subsample.grid,
               xgb.colsample.grid = xgb.colsample.grid,
               cb.iterations.grid = cb.iterations.grid,
               cb.depth.grid = cb.depth.grid,
               cb.learning_rate.grid = cb.learning_rate.grid,
               cb.l2_leaf_reg.grid = cb.l2_leaf_reg.grid,
               torch.hidden_units.grid = torch.hidden_units.grid,
               torch.learning_rate.grid = torch.learning_rate.grid,
               torch.dropout.grid = torch.dropout.grid,
               torch.epochs.grid = torch.epochs.grid)
}
