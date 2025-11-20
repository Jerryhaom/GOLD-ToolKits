#' GOLD-R Biological Age Estimation
#'
#' Calculate biological age using the GOLD-R method based on residual prediction
#' from biomarker data.
#'
#' @param d4 A data frame containing the input data
#' @param var Character vector of variable names to use
#' @param feature_selection Logical indicating whether to perform feature selection
#' @param selection_method Method for feature selection ("lasso", "ridge", "elasticnet")
#' @param alpha Elasticnet mixing parameter (0 for ridge, 1 for lasso)
#' @param nfolds Number of folds for cross-validation
#' @param family Family for glmnet ("gaussian", "cox", "binomial")
#'
#' @return A list containing:
#' \item{GOLDR-Bioage}{Estimated biological ages}
#' \item{residual}{Model residuals}
#' \item{coef}{Model coefficients}
#' \item{age}{Chronological age}
#' \item{biomarker_vars}{Biomarker variables used}
#'
#' @examples
#' \dontrun{
#' Demo example on simulated datasets
#'data1=generate_data_sim()
#' var_cols <- setdiff(colnames(d4), c("id", "time","status"))
#' var=var_cols
#' result <- goldr_bioage(
#'  d4 = d4,
#'  var = var,
#'  feature_selection = TRUE,
#'  selection_method = "lasso",
#'  alpha = 1,
#'  nfolds = 10,
#'  family="gaussian"
#')
#'names(result)
#'
#' @export
goldr_bioage <- function(d4, var, feature_selection = FALSE,
                         selection_method = "lasso", alpha = 1,
                         nfolds = 10, family = "gaussian") {

  # Validate inputs
  if (!"age" %in% var) {
    stop("'age' must be included in the variable list")
  }

  if (!all(c("time", "status") %in% names(d4))) {
    stop("Data must contain 'time' and 'status' columns")
  }

  # Prepare data
  data <- data.frame(d4[, var, drop = FALSE], check.names = FALSE)
  data <- stats::na.omit(data)

  # Dynamic identification of biomarker variables
  biomarker_vars <- setdiff(var, c("age", "time", "status", "id"))

  if (length(biomarker_vars) == 0) {
    stop("No biomarker variables found. Please include biomarker columns in 'var'")
  }

  # Step 1: calculate initial residuals from predicting Age-related mortality risk
  x_data <- as.matrix(data[, biomarker_vars, drop = FALSE])
  shape <- 0.08
  y_data <- shape * data[, "age"] + log(0.00002)

  # Perform cross-validated LASSO
  set.seed(123)  # for reproducibility
  cv_fit <- glmnet::cv.glmnet(
    x = x_data,
    y = y_data,
    family = family,
    alpha = alpha,
    nfolds = nfolds
  )

  # Predict and calculate residuals
  hhat <- predict(cv_fit, x_data)
  hhat <- as.matrix(hhat)
  dd <- data.frame(hhat = hhat[, 1], age = data[, "age"])
  l <- stats::lm(hhat ~ age, data = dd)
  a1 <- stats::coef(l)[1]
  b1 <- stats::coef(l)[2]
  s <- summary(l)
  resi <- hhat - b1 * dd$age

  # Step 2: predict residuals using biomarkers
  cv_fit2 <- glmnet::cv.glmnet(
    x = x_data,
    y = resi,
    family = family,
    alpha = alpha,
    nfolds = nfolds
  )

  yhat <- predict(cv_fit2, x_data)
  yhat <- as.matrix(yhat)

  # Scale residuals - 修复：直接使用mean()而不是stats::mean()
  resi1 <- yhat[, 1] * stats::sd(resi) / stats::sd(yhat[, 1])
  gamma <- stats::sd(resi) / stats::sd(yhat[, 1])
  beta0 <- -mean(resi1) + mean(resi)  # 修复：直接使用mean()
  resi1 <- resi1 + beta0

  # Calculate biological age adjustment
  dage <- (resi1 - s$coefficients[1, 1]) / s$coefficients[2, 1]

  # Extract coefficients
  coef_matrix <- coef(cv_fit2)
  coef2 <- as.numeric(coef_matrix)
  names(coef2) <- rownames(coef_matrix)
  coef2 <- coef2[coef2 != 0]

  # Adjust coefficients
  coef3 <- coef2 * gamma
  if ("(Intercept)" %in% names(coef3)) {
    coef3["(Intercept)"] <- coef3["(Intercept)"] - s$coefficients[1, 1] + beta0
  }
  coef3 <- coef3 / s$coefficients[2, 1]

  # Prepare results
  result <- list(
    residual = dage,
    coef = coef3,
    `GOLDR-Bioage` = data[, "age"] + dage,
    age = data[, "age"],
    biomarker_vars = biomarker_vars
  )

  return(result)
}
