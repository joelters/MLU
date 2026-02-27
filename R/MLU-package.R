#' MLU: Machine Learning U-Statistics
#'
#' The MLU package provides functions for estimation of conditional U-statistics 
#' with machine learning for pairwise data structures. It includes tools 
#' for model estimation, hyperparameter tuning, and fitted value 
#' estimation for pairwise data relationships.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{mlumodest}}: Model estimation for pairwise data
#'   \item \code{\link{mlutuning}}: Hyperparameter tuning across multiple algorithms
#'   \item \code{\link{mluFVest}}: Fitted value estimation with treatment assignment
#'   \item \code{\link{mluFVestab}}: Fitted values for all treatment combinations
#' }
#'
#' @section Pairwise Data:
#' Pairwise data consists of observations on pairs of units. These functions create
#' all possible pairs (i,j) where i < j ("triangle" mode) or all combinations 
#' between two groups ("square" mode), and model the relationship using machine 
#' learning methods from the ML package.
#'
#' @keywords internal
"_PACKAGE"
