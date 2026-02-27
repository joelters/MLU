# MLU: Machine Learning U-Statistics

## Overview

MLU provides functions for estimation of conditional U-statistics with machine learning for dyadic (pairwise) data structures. The package implements methods for analyzing data where observations are paired, allowing researchers to model relationships between pairs of units using state-of-the-art machine learning techniques.

## Installation

You can install the development version of MLU from this directory:

```r
# Install from local source
install.packages("path/to/MLU", repos = NULL, type = "source")

# Or using devtools
devtools::install_local("path/to/MLU")
```

## Main Functions

### Model Estimation

- `dyadmodest()`: Estimates machine learning models for dyadic data by creating all possible pairs (i,j) where i < j and fitting a model using combined features from both observations.

### Cross-Validation and Tuning

- `dyadcv()`: Performs cross-validation to select the best machine learning algorithm based on RMSE.
- `dyadtuning()`: Comprehensive hyperparameter tuning across multiple ML algorithms including:
  - Random Forest (RF)
  - Conditional Inference Forest (CIF)
  - XGBoost (XGB)
  - CatBoost (CB)
  - Lasso and Ridge regression
  - Logistic Lasso
  - Neural networks (torch)

### Fitted Value Estimation

- `dyadFVest()`: Estimates fitted values for dyadic models with two modes:
  - "triangle": within-group pairs (i < j from same group)
  - "square": between-group pairs (all i,j combinations from two groups)

- `dyadFVestab()`: Estimates fitted values for all treatment assignment combinations (1,1), (1,0), (0,1), (0,0), useful for causal inference with dyadic data.

## Basic Usage

```r
library(MLU)

# Example data
n <- 100
X <- matrix(rnorm(n * 5), n, 5)
D <- rbinom(n, 1, 0.5)
Y <- rnorm(n)

# Define a dyadic function (e.g., absolute difference)
f <- function(y1, y2) abs(y1 - y2)

# Fit a dyadic model
model <- dyadmodest(D, X, Y, f, ML = c("Lasso", "RF"))

# Cross-validation to select best learner
cv_results <- dyadcv(X, Y, f, ML = c("Lasso", "Ridge", "RF"), Kcv = 5)
print(cv_results$mlbest)  # Best performing algorithm
print(cv_results$rmse)     # RMSE of best algorithm

# Hyperparameter tuning for specific algorithms
tuning_results <- dyadtuning(X, Y, f, ML = c("RF", "XGB"), Kcv = 5)

# Estimate fitted values (within-group pairs)
fv <- dyadFVest(model, D, X, Y, D, X, Y, f, shape = "triangle")

# Estimate fitted values for all treatment combinations
fv_all <- dyadFVestab(model, X, Y, X, Y, f, shape = "triangle")
print(fv_all$fv11)  # Both treated
print(fv_all$fv00)  # Both untreated
```

## Dependencies

This package requires the `ML` package, which provides the underlying machine learning infrastructure.

## Source

The dyadic functions in this package are adapted from the [ineqewm](https://github.com/joelters/ineqewm) package by Joel Terrier, which implements methods for inequality-efficient welfare maximization.

## License

MIT

## References

For methodological details on dyadic data analysis and U-statistics with machine learning, see the documentation and references in the original ineqewm package.
