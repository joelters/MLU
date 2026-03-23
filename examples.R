################################################################################
# MLU Package Examples
# This file demonstrates basic usage of the MLU package functions
################################################################################

# Load the package (after installation)
library(MLU)

################################################################################
# Example 1: Basic pairwise model estimation
################################################################################

# Generate example data
set.seed(123)
n <- 100
X <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n)
)
D <- rbinom(n, 1, 0.5)  # Treatment assignment
Y <- 2 * D + X$x1 + 0.5 * X$x2 + rnorm(n, sd = 0.5)  # Outcome

# Define a pairwise function (e.g., absolute difference)
f_abs_diff <- function(y1, y2) abs(y1 - y2)

# Fit a pairwise model
model <- mlumodest(D, X, Y, f = f_abs_diff, ML = c("OLS"), polynomial.OLS = 2)
print("Model fitted successfully!")

################################################################################
# Example 2: Hyperparameter tuning
################################################################################

# Tune hyperparameters for specific algorithms
# Note: This can be computationally intensive
tuning_results <- mlutuning(
  X = X,
  Y = Y,
  f = f_abs_diff,
  ML = c("RF"),
  Kcv = 3,
  rf.cf.ntree.grid = c(100, 200),
  rf.depth.grid = c(2, 4),
  mtry.grid = c(1, 2),
  subsample = NULL  # Use 50 observations for faster computation
)

cat("\nHyperparameter tuning completed!\n")

################################################################################
# Example 4: Fitted value estimation (triangle mode)
################################################################################

# Estimate fitted values for within-group pairs
# Create a new dataset for prediction
n_new <- 20
X_new <- data.frame(
  x1 = rnorm(n_new),
  x2 = rnorm(n_new),
  x3 = rnorm(n_new)
)
D_new <- rbinom(n_new, 1, 0.5)
Y_new <- 2 * D_new + X_new$x1 + 0.5 * X_new$x2 + rnorm(n_new, sd = 0.5)

fv <- mluFVest(
  model = model,
  Di = D,
  Xi = X,
  Yi = Y,
  Dnewi = D_new,
  Xnewi = X_new,
  Ynewi = Y_new,
  f = f_abs_diff,
  ML = "OLS"
)

cat("\nNumber of fitted values (triangle):", length(fv), "\n")
cat("Expected number of pairs:", n_new * (n_new - 1) / 2, "\n")

################################################################################
# Example 5: Fitted values for all treatment combinations
################################################################################

# Estimate fitted values for all treatment assignment combinations
# Useful for causal inference
fv_all <- mluFVestab(
  model = model,
  Xi = X,
  Yi = Y,
  Xnewi = X_new[1, ],
  Ynewi = Y_new[1],
  Xj = X,
  Yj = Y,
  Xnewj = X_new[2, ],
  Ynewj = Y_new[2],
  f = f_abs_diff,
  ML = "OLS",
  polynomial.OLS = 2,
  shape = "square"
)

cat("\nFitted values for all treatment combinations:\n")
cat("  fv11 (both treated):", length(fv_all$fv11), "values\n")
cat("  fv10 (i treated, j untreated):", length(fv_all$fv10), "values\n")
cat("  fv01 (i untreated, j treated):", length(fv_all$fv01), "values\n")
cat("  fv00 (both untreated):", length(fv_all$fv00), "values\n")

################################################################################
# Example 6: Using different pairwise functions
################################################################################

# Define alternative pairwise functions
f_min <- function(y1, y2) pmin(y1, y2)
f_max <- function(y1, y2) pmax(y1, y2)
f_mean <- function(y1, y2) (y1 + y2) / 2
f_gini <- function(y1, y2) 0.5 * (y1 + y2 - abs(y1 - y2))  # Gini-like

# Fit models with different pairwise functions
model_min <- mlumodest(D, X, Y, f = f_min, ML = "Lasso")
model_max <- mlumodest(D, X, Y, f = f_max, ML = "Lasso")
model_mean <- mlumodest(D, X, Y, f = f_mean, ML = "Lasso")
model_gini <- mlumodest(D, X, Y, f = f_gini, ML = "Lasso")

cat("\nAll pairwise function models fitted successfully!\n")

################################################################################
# Example 7: Square mode (between-group pairs)
################################################################################

# Split data into two groups
n1 <- 15
n2 <- 15
X_group1 <- X[1:n1, ]
X_group2 <- X[(n1+1):(n1+n2), ]
D_group1 <- D[1:n1]
D_group2 <- D[(n1+1):(n1+n2)]
Y_group1 <- Y[1:n1]
Y_group2 <- Y[(n1+1):(n1+n2)]

# Estimate fitted values for between-group pairs
fv_square <- mluFVest(
  model = model,
  Di = D_group1,
  Xi = X_group1,
  Yi = Y_group1,
  Dnewi = D_group1,
  Xnewi = X_group1,
  Ynewi = Y_group1,
  Dj = D_group2,
  Xj = X_group2,
  Yj = Y_group2,
  Dnewj = D_group2,
  Xnewj = X_group2,
  Ynewj = Y_group2,
  f = f_abs_diff,
  ML = "Lasso",
  shape = "square"
)

cat("\nNumber of fitted values (square):", length(fv_square), "\n")
cat("Expected number of pairs:", n1 * n2, "\n")

cat("\n################################################################################\n")
cat("All examples completed successfully!\n")
cat("################################################################################\n")

################################################################################
# Example 8: Explicit pairwise reference example
################################################################################

# Separate example with 100 observations to illustrate the cached pairwise
# reference data created by mlumodest().
set.seed(321)
n_cache <- 100
D_cache <- rbinom(n_cache, 1, 0.5)
X_cache <- data.frame(
  x1 = rnorm(n_cache),
  x2 = rnorm(n_cache)
)
Y_cache <- 1 + 2 * D_cache + X_cache$x1 - 0.5 * X_cache$x2 + rnorm(n_cache, sd = 0.25)

# Pairwise outcome function used for illustration
f_cache <- function(y1, y2) abs(y1 - y2)

# Fit the model; mlumodest now stores the fitted learner plus cached
# pairwise training reference data used by mluFVest / mluFVestab.
model_cache <- mlumodest(
  D_cache,
  X_cache,
  Y_cache,
  f = f_cache,
  ML = "OLS",
  polynomial.OLS = 2,
  verbose = FALSE
)

cat("\nCached model object contents:\n")
print(names(model_cache))

cat("\nCached pairwise training reference dimensions:\n")
print(dim(model_cache$Xref))

cat("\nFirst 6 rows of cached pairwise training reference Xref:\n")
print(utils::head(model_cache$Xref))

cat("\nFirst 6 cached pairwise outcomes Yref:\n")
print(utils::head(model_cache$Yref))

cat("\nFor 100 observations there are 100 * 99 / 2 =", n_cache * (n_cache - 1) / 2,
    "pairwise training rows cached in Xref.\n")
cat("Each Xref row stores [D_i + D_j, X_i + X_j, |D_i - D_j|, |X_i - X_j|].\n")

# One new pair for prediction in square mode
X_cache_new_i <- data.frame(x1 = 0.25, x2 = -0.50)
Y_cache_new_i <- 1.8
X_cache_new_j <- data.frame(x1 = -1.10, x2 = 0.75)
Y_cache_new_j <- 0.2

fv_cache <- mluFVestab(
  model = model_cache,
  Xi = X_cache,
  Yi = Y_cache,
  Xnewi = X_cache_new_i,
  Ynewi = Y_cache_new_i,
  Xj = X_cache,
  Yj = Y_cache,
  Xnewj = X_cache_new_j,
  Ynewj = Y_cache_new_j,
  f = f_cache,
  ML = "OLS",
  polynomial.OLS = 2,
  shape = "square"
)

cat("\nOne-pair square predictions from cached reference data:\n")
print(fv_cache)