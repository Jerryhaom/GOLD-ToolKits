#' Generate Simulated Survival Data with Biomarkers
#'
#' @param n Sample size
#' @param p Number of biomarkers
#' @param gamma_val Shape parameter for Gompertz distribution
#' @param lb_val Scale parameter for Gompertz distribution
#' @param num_signal Number of true signal biomarkers
#' @param sigma_age Standard deviation for age noise
#' @param beta_age Effect size for age
#' @param beta_extra Effect size for extra signal
#' @param maxt Maximum follow-up time
#' @param age_min Minimum age
#' @param age_max Maximum age
#' @param beta_min_val Minimum effect size for biomarkers
#' @param beta_max_val Maximum effect size for biomarkers
#' @return A data frame with simulated survival data
#' @importFrom MASS mvrnorm
#' @importFrom simsurv simsurv
#' @importFrom stats toeplitz runif rnorm
#' @export
generate_data_sim <- function(
    n = 1000,        # 减小默认值便于测试
    p = 50,         # 减小默认值便于测试
    gamma_val = 0.1,
    dt = "gompertz",
    lb_val = 1e-06,
    num_signal = 10,
    sigma_age = 5,
    beta_age = 0.1,
    beta_extra = 1,
    maxt = 8,
    age_min = 30,
    age_max = 70,
    beta_min_val = 0.05,
    beta_max_val = 0.15
) {

  # 设置随机种子以确保可重复性
  set.seed(123)

  # 生成 biomarker 数据矩阵 - 修复参数名
  mu <- rep(0, p)
  rho <- 0.3
  sigma_auto <- stats::toeplitz(rho^(0:(p-1)))
  diag(sigma_auto) <- 1

  # 修复：使用正确的参数名 mu 和 Sigma
  data_matrix <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma_auto)
  df <- as.data.frame(scale(data_matrix))
  colnames(df) <- paste0("X", 1:p)

  # 随机选择信号 biomarker
  signal_idx <- sample(1:p, num_signal)

  # 生成 age 对应的 beta
  true_betas <- rep(0, p)
  true_betas[signal_idx] <- stats::runif(num_signal, min = beta_min_val, max = beta_max_val)

  # 生成 Age = β·X + ε
  y_values <- as.matrix(df) %*% true_betas

  # 归一化 age 到用户指定区间
  age_values <- (y_values - min(y_values)) / (max(y_values) - min(y_values)) *
    (age_max - age_min) + age_min + stats::rnorm(n, 0, sigma_age)

  # 准备协变量数据
  covdat <- data.frame(age = age_values)

  # 生成额外的信号
  signal_idx2 <- sample(1:p, num_signal)
  true_betas2 <- rep(0, p)
  true_betas2[signal_idx2] <- stats::runif(num_signal, min = beta_min_val, max = beta_max_val)
  z_values <- as.matrix(df) %*% true_betas2
  z_values <- (z_values - min(z_values)) / (max(z_values) - min(z_values))
  covdat$extra <- z_values

  # 设置效应系数
  betas_final <- c(age = beta_age, extra = beta_extra)

  # 生成生存数据
  surv_data <- simsurv::simsurv(
    dist = dt,
    lambdas = lb_val,
    gammas = gamma_val,
    x = covdat,
    betas = betas_final,
    maxt = maxt
  )

  # 输出完整数据
  result <- data.frame(
    id = 1:n,
    df,
    time = surv_data$eventtime,
    status = surv_data$status,
    age = age_values
  )

  return(result)
}
