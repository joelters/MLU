#' @keywords internal
"_PACKAGE"

#' MLU: Machine Learning U-Statistics
#'
#' The MLU package provides functions for estimation of conditional U-statistics 
#' with machine learning for dyadic (pairwise) data structures. It includes tools 
#' for model estimation, cross-validation, hyperparameter tuning, and fitted value 
#' estimation for pairwise data relationships.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{dyadmodest}}: Model estimation for dyadic data
#'   \item \code{\link{dyadcv}}: Cross-validation for algorithm selection
#'   \item \code{\link{dyadtuning}}: Hyperparameter tuning across multiple algorithms
#'   \item \code{\link{dyadFVest}}: Fitted value estimation with treatment assignment
#'   \item \code{\link{dyadFVestab}}: Fitted values for all treatment combinations
#' }
#'
#' @section Dyadic Data:
#' Dyadic data consists of observations on pairs of units. These functions create
#' all possible pairs (i,j) where i < j ("triangle" mode) or all combinations 
#' between two groups ("square" mode), and model the relationship using machine 
#' learning methods from the ML package.
#'
#' @docType package
#' @name MLU-package
NULL
