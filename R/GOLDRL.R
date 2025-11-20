
#' Calculate biological age using Gompertz model with feature selection
#'
#' @description
#' Computes biological age a Gompertz regression model
#' Cox-LASSO feature selection for biomarkers
#'
#' @param d4 Dataframe containing survival information (time, status)
#' @param var Character vector of variable names (must include 'age')
#' @param feature_selection Logical indicating whether to perform feature selection (default: FALSE)
#' @param selection_method Method for feature selection: "lasso", "elasticnet", or "none"
#' @param alpha Elasticnet mixing parameter (0 = ridge, 1 = lasso)
#' @param nfolds Number of folds for cross-validation
#' @param family Survival family for glmnet (default: "cox")
#'
#' @return A list containing:
#' \itemize{
#'   \item residual: GOLD-RL BioAge Residuals
#'   \item Coef: Coefficients of Biomarkers
#'   \item GOLDR-Bioage: GOLD-RL BioAge
#'   \item risk: moratlity hazard from Gompertz regression model
#' }
#'
#' @importFrom stats as.formula na.omit predict
#' @importFrom stats lm coef rnorm
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs theme_minimal
#' @export
#'
#' @examples
#' # Load dataset
#' data(NHANES4)
#'
#' # Define variables
#' var <- c("age", "albumin", "alp", "creat", "glucose_mmol", "lymph", "mcv", "rdw", "wbc", "ggt")
#'
#' # Without feature selection
#' result1 <- goldrl_bioage(NHANES4, var)
#'
#' # With LASSO feature selection
#' result2 <- goldrl_bioage(NHANES4, var,
#'                       feature_selection = TRUE,
#'                       selection_method = "lasso")
#'
#' # With Elasticnet feature selection
#' result3 <- goldrl_bioage(NHANES4, var,
#'                       feature_selection = TRUE,
#'                       selection_method = "elasticnet",
#'                       alpha = 0.5)
#'
##GOLD-RL
goldrl_bioage <- function(d4, var, feature_selection = TRUE,
                          selection_method = "lasso", alpha = 1,
                          nfolds = 5, family = "cox") {
  # Validate inputs
  if (!"age" %in% var) {
    stop("'age' must be included in the variable list")
  }
  if (!all(c("time", "status") %in% names(d4))) {
    stop("Data must contain 'time' and 'status' columns")
  }
  # Prepare data
  data <- data.frame(d4[, c("time", "status", var)], check.names = FALSE)
  data <- data[data$time > 0, ]
  data <- na.omit(data)
  colnames(data) <- make.names(colnames(data))
  # Store original variables
  original_vars <- var
  biomarker_vars <- setdiff(var, "age")
  # Feature selection if requested
  selected_vars <- NULL
  cv_plot <- NULL
  if (feature_selection && length(biomarker_vars) > 1) {
    # Prepare data for glmnet
    x_data <- as.matrix(data[, make.names(biomarker_vars)])
    y_data <- survival::Surv(data$time, data$status)
    # Set alpha based on selection method
    if (selection_method == "lasso") {
      alpha_val <- 1
    } else if (selection_method == "elasticnet") {
      alpha_val <- alpha
    } else {
      alpha_val <- 1  # default to lasso
    }
    # Perform cross-validated LASSO
    set.seed(123)  # for reproducibility
    cv_fit <- glmnet::cv.glmnet(x = x_data, y = y_data,
                                family = family,
                                alpha = alpha_val,
                                nfolds = nfolds)
    # Get selected variables (non-zero coefficients at lambda.min)
    coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
    selected_biomarkers <- biomarker_vars[coefs[-1, 1] != 0]  # exclude intercept

    # Always include Age and selected biomarkers
    selected_vars <- c("age", make.names(selected_biomarkers))
    selected_vars <- selected_vars[selected_vars %in% colnames(data)]  # 安全过滤
    message(paste("Selected", length(selected_biomarkers),
                  "biomarkers out of", length(biomarker_vars), "using", selection_method))
  } else {
    # No feature selection or not enough biomarkers
    selected_vars <- make.names(var)
    if (feature_selection && length(biomarker_vars) <= 1) {
      warning("Not enough biomarkers for feature selection. Using all variables.")
    }
  }
  data1 <- data
  # Model: Age + (selected) biomarkers
  formula_str <- paste(
    "survival::Surv(time, status) ~",
    paste(selected_vars, collapse = "+")
  )
  fitg2 <- flexsurv::flexsurvreg(
    formula = as.formula(formula_str),
    data = data1, dist = "gompertz"
  )
  coef2 <- fitg2$res[, 1]
  # Calculate mortality risk
  pr <- predict(fitg2, type = "hazard", times = c(0))
  risk <- pr$.pred_hazard

  hhat=log(risk)

  # Step 1: calculate initial residuals for predicting mortality risk
  dd = data.frame(hhat = log(risk),
                  age = data[,"age"])
  l <- lm(hhat ~ age, data = dd)
  a1=coef(l)[1]
  b1=coef(l)[2]
  s=summary(l)
  resi=hhat-b1*dd$age
  # Step 1: predict residuals using rest biomarkers
  cv_fit2 <- glmnet::cv.glmnet(x = x_data, y = resi, family = "gaussian",alpha = alpha_val,nfolds = nfolds)
  yhat = predict(cv_fit2, x_data)
  yhat <- as.matrix(yhat)
  resi1 = yhat[,1] * sd(resi) / sd(yhat[,1])
  gamma = sd(resi) / sd(yhat[,1])
  beta0 = -mean(resi1) + mean(resi)
  resi1 = resi1 + beta0
  dage = (resi1 - s$coefficients[1,1]) / s$coefficients[2,1]
  coef2 = coef(cv_fit2)[,1]
  coef2 = coef2[coef2 != 0]
  coef3 = coef2 * gamma
  coef3[1] = coef3[1] - s$coefficients[1,1] + beta0
  coef3 = coef3 / s$coefficients[2,1]
  # Prepare results
  result <- list(
    d4,
    residual = dage,
    coef = coef3,
    `GOLDR-Bioage` = data[,"age"] + dage,
    risk=risk
  )
  return(result)
}
