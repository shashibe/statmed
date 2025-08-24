#------------------------
rm(list = ls())
set.seed(1234)

# ---- packages ----
library(MASS)
library(logistf)
library(geepack)
library(detectseparation)
library(parallel)
library(pbapply)

# ---- directories ----
data_dir <- "data"       # 入力: 真値のCSVを置くフォルダ
out_dir  <- "results"    # 出力: シミュレーション結果を保存するフォルダ
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- argument ----
N_sim    <- 1000
n        <- c(2500, 500)
rho      <- c(0.5, 0.25, 0)
event    <- c(0.005, 0.01, 0.03, 0.05, 0.1, 0.2)
exposure <- 0.1
beta0    <- c(
  -12.9, -12.2, -10.85, -10.25, -9.1, -7.6,
  -12.7, -11.5, -10.5, -9.6, -8.5, -7.3,
  -11.75, -11.2, -10.05, -9.3, -8.25, -7.1,
  -12.24, -11.6, -10.9, -10, -8.9, -7.55,
  -11.96, -11.5, -10.3, -9.7, -8.7, -7.2,
  -11.8, -11.3, -10.05, -9.35, -8.3, -7.03
) 
alpha0   <- c(
  -8.9, -8.9, -8.9, -8.9, -8.9, -8.9,
  -8.4, -8.4, -8.4, -8.4, -8.4, -8.4,
  -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,
  -8.55, -8.55, -8.55, -8.55, -8.55, -8.55,
  -8.2, -8.2, -8.2, -8.2, -8.2, -8.2,
  -8.2, -8.2, -8.2, -8.2, -8.2, -8.2
)

# ---- scenario summary ----
rows <- list()
for (N in n) {
  for (cor in rho) {
    for (ev in event) {
      rows[[length(rows) + 1]] <- data.frame(
        n = N,
        event = ev,
        exposure = exposure,
        rho = cor
      )
    }
  }
}
sim_summary <- do.call(rbind, rows)
rownames(sim_summary) <- NULL
sim_summary <- cbind(sim_summary, beta0 = beta0, alpha0 = alpha0, N_sim = N_sim)

# ---- simulation function ----
run_simulation <- function(row_index) {
  scenario <- sim_summary[row_index, ]
  
  # ---- read true value CSV ----
  true_value <- read.csv(file.path(data_dir, paste0(
    "sim-true-value-rho=", scenario$rho,
    "-event=", scenario$event,
    "-exposure=", scenario$exposure, "-20230610.csv"
  )))
  
  mu1       <- true_value$mu_1
  mu0       <- true_value$mu_0
  diff      <- true_value$diff
  log_ratio <- true_value$log_ratio
  
  #----Estimation of causal effects
  mrg_rsk_cal <- function(est, cov_1, cov_0){
    if(length(est) <= 2){
      xbeta_1 <- est[1] + est[2]
      xbeta_0 <- est[1] 
      rsk_1 <- mean(exp(xbeta_1) / (1 + exp(xbeta_1)))
      rsk_0 <- mean(exp(xbeta_0) / (1 + exp(xbeta_0)))
      rsk_diff <- rsk_1 - rsk_0
      rsk_log_ratio <- log(rsk_1 / rsk_0)
      rsk_ratio     <- rsk_1 / rsk_0
    }else{
      xbeta_1 <- cov_1 %*% est
      xbeta_0 <- cov_0 %*% est
      rsk_1 <- mean(1 / (1 + exp(-xbeta_1)))
      rsk_0 <- mean(1 / (1 + exp(-xbeta_0)))
      rsk_diff <- rsk_1 - rsk_0
      rsk_log_ratio <- log(rsk_1 / rsk_0)
      rsk_ratio     <- rsk_1 / rsk_0
    }
    df <- data.frame(rsk_1 = rsk_1, rsk_0 = rsk_0, rsk_diff = rsk_diff,
                     rsk_log_ratio = rsk_log_ratio, rsk_ratio = rsk_ratio)
    return(df)
  }
  #----Calcuration of coverage_count
  coverage_count <- function(est, cov_1, cov_0, covariance, risk_1, risk_0, risk_diff, risk_log_ratio, n){
    if(is.na(risk_1) == T || is.na(risk_0) == T){
      count_0         <- 0
      count_1         <- 0
      count_diff      <- 0
      count_log_ratio <- 0
      se_1            <- NA
      se_0            <- NA
      se_diff         <- NA
      se_log_ratio    <- NA
      df <- data.frame(count_1 = count_1, count_0 = count_0, count_diff = count_diff, count_log_ratio = count_log_ratio,
                       se_1 = se_1, se_0 = se_0, se_diff = se_diff, se_log_ratio = se_log_ratio)
    }else{
      # asymptotic standard error
      xbeta_1 <- cov_1 %*% est
      xbeta_0 <- cov_0 %*% est
      A_1 <- exp(xbeta_1) / (1 + exp(xbeta_1)) ^ 2
      A_0 <- exp(xbeta_0) / (1 + exp(xbeta_0)) ^ 2
      
      D_1         <- as.matrix((as.vector(A_1) %*% cov_1) / n)
      D_0         <- as.matrix((as.vector(A_0) %*% cov_0) / n)
      D_diff      <- as.matrix((as.vector(A_1) %*% cov_1 - as.vector(A_0) %*% cov_0) / n)
      D_log_ratio <- as.matrix((as.vector(A_1) %*% cov_1) / (risk_1 * n ) - (as.vector(A_0) %*% cov_0) / (risk_0 * n ))
      
      se_1         <- sqrt(D_1 %*% covariance %*% t(D_1))
      se_0         <- sqrt(D_0 %*% covariance %*% t(D_0))
      se_diff      <- sqrt(D_diff %*% covariance %*% t(D_diff))
      se_log_ratio <- sqrt(D_log_ratio %*% covariance %*% t(D_log_ratio))
      
      upper_1         <- risk_1 + 1.96 * se_1
      lower_1         <- risk_1 - 1.96 * se_1
      upper_0         <- risk_0 + 1.96 * se_0
      lower_0         <- risk_0 - 1.96 * se_0
      upper_diff      <- risk_diff + 1.96 * se_diff
      lower_diff      <- risk_diff - 1.96 * se_diff
      upper_log_ratio <- risk_log_ratio + 1.96 * se_log_ratio
      lower_log_ratio <- risk_log_ratio - 1.96 * se_log_ratio
      
      count_1         <- ifelse(lower_1 < risk_1 && risk_1 < upper_1, 1, 0)
      count_0         <- ifelse(lower_0 < risk_0 && risk_0 < upper_0, 1, 0)
      count_diff      <- ifelse(lower_diff < risk_diff && risk_diff < upper_diff, 1, 0)
      count_log_ratio <- ifelse(lower_log_ratio < risk_log_ratio && risk_log_ratio < upper_log_ratio, 1, 0)
      df <- data.frame(count_1 = count_1, count_0 = count_0, count_diff = count_diff, count_log_ratio = count_log_ratio,
                       se_1 = se_1, se_0 = se_0, se_diff = se_diff, se_log_ratio = se_log_ratio,
                       upper_1 = upper_1, upper_0 = upper_0, upper_diff = upper_diff, upper_log_ratio = upper_log_ratio,
                       lower_1 = lower_1, lower_0 = lower_0, lower_diff = lower_diff, lower_log_ratio = lower_log_ratio)
    }
    
    return(df)
  }
  
  #---- result data
  res_rsk_1 <- data.frame(n = scenario$n,
                          measure = "mu1",
                          Itr = 1 : scenario$N_sim,
                          sepa = 0,
                          ML = 0,
                          no_adj_ML = 0,
                          FML = 0, 
                          FLIC = 0,
                          FLAC = 0,
                          ps_adj = 0,
                          IPW = 0)
  res_rsk_0 <- data.frame(n = scenario$n,
                          measure = "mu0",
                          Itr = 1 : scenario$N_sim,
                          sepa = 0,
                          ML = 0,
                          no_adj_ML = 0,
                          FML = 0, 
                          FLIC = 0,
                          FLAC = 0,
                          ps_adj = 0,
                          IPW = 0)
  res_rsk_diff <- data.frame(n = scenario$n,
                             measure = "diff",
                             Itr = 1 : scenario$N_sim,
                             sepa = 0,
                             ML = 0,
                             no_adj_ML = 0,
                             FML = 0, 
                             FLIC = 0,
                             FLAC = 0,
                             ps_adj = 0,
                             IPW = 0)
  res_rsk_log_ratio <- data.frame(n = scenario$n,
                                  measure = "log_ratio",
                                  Itr = 1 : scenario$N_sim,
                                  sepa = 0,
                                  ML = 0,
                                  no_adj_ML = 0,
                                  FML = 0, 
                                  FLIC = 0,
                                  FLAC = 0,
                                  ps_adj = 0,
                                  IPW = 0)
  coverage_count_rsk_1 <- data.frame(n = scenario$n,
                                     Itr = 1 : scenario$N_sim,
                                     sepa = 0,
                                     ML = 0,
                                     no_adj_ML = 0,
                                     FML = 0, 
                                     FLIC = 0,
                                     FLAC = 0,
                                     ps_adj = 0,
                                     IPW = 0)
  coverage_count_rsk_0 <- data.frame(n = scenario$n,
                                     Itr = 1 : scenario$N_sim,
                                     sepa = 0,
                                     ML = 0,
                                     no_adj_ML = 0,
                                     FML = 0, 
                                     FLIC = 0,
                                     FLAC = 0,
                                     ps_adj = 0,
                                     IPW = 0)
  coverage_count_rsk_diff <- data.frame(n = scenario$n,
                                        Itr = 1 : scenario$N_sim,
                                        sepa = 0,
                                        ML = 0,
                                        no_adj_ML = 0,
                                        FML = 0, 
                                        FLIC = 0,
                                        FLAC = 0,
                                        ps_adj = 0,
                                        IPW = 0)
  coverage_count_rsk_log_ratio <- data.frame(n = scenario$n,
                                             Itr = 1 : scenario$N_sim,
                                             sepa = 0,
                                             ML = 0,
                                             no_adj_ML = 0,
                                             FML = 0, 
                                             FLIC = 0,
                                             FLAC = 0,
                                             ps_adj = 0,
                                             IPW = 0)
  mese_count_rsk_1 <- data.frame(n = scenario$n,
                                 Itr = 1 : scenario$N_sim,
                                 sepa = 0,
                                 ML = 0,
                                 no_adj_ML = 0,
                                 FML = 0, 
                                 FLIC = 0,
                                 FLAC = 0,
                                 ps_adj = 0,
                                 IPW = 0)
  mese_count_rsk_0 <- data.frame(n = scenario$n,
                                 Itr = 1 : scenario$N_sim,
                                 sepa = 0,
                                 ML = 0,
                                 no_adj_ML = 0,
                                 FML = 0, 
                                 FLIC = 0,
                                 FLAC = 0,
                                 ps_adj = 0,
                                 IPW = 0)
  mese_count_rsk_diff <- data.frame(n = scenario$n,
                                    Itr = 1 : scenario$N_sim,
                                    sepa = 0,
                                    ML = 0,
                                    no_adj_ML = 0,
                                    FML = 0, 
                                    FLIC = 0,
                                    FLAC = 0,
                                    ps_adj = 0,
                                    IPW = 0)
  mese_count_rsk_log_ratio <- data.frame(n = scenario$n,
                                         Itr = 1 : scenario$N_sim,
                                         sepa = 0,
                                         ML = 0,
                                         no_adj_ML = 0,
                                         FML = 0, 
                                         FLIC = 0,
                                         FLAC = 0,
                                         ps_adj = 0,
                                         IPW = 0)
  for (i in 1 : scenario$N_sim) {
    set.seed(1234 + i)
    
    # multivariate normal distribution
    beta     <- c( rep(1, 10))
    alpha    <- c( rep(1, 10))
    betaz    <- 0.8
    mu         <- c(rep(0, length(beta)))
    sigma      <- c(rep(1, length(beta)))
    Rho        <- scenario$rho ^ abs(outer(1:10, 1:10, "-"))
    sigma_diag <- diag(sigma)
    Sigma      <- sigma_diag %*% Rho %*% sigma_diag
    multi_X    <- mvrnorm(scenario$n, mu, Sigma)
    
    # generate covariate
    X <- cbind(
      x1  <- as.numeric(0 <= multi_X[, 1]),
      x2  <- as.numeric(0 <= multi_X[, 2]),
      x3  <- as.numeric(0 <= multi_X[, 3]),
      x4  <- as.numeric(0 <= multi_X[, 4]),
      x5  <- as.numeric(0 <= multi_X[, 5]),
      x6  <- as.numeric(0 <= multi_X[, 6]),
      x7  <- as.numeric(0 <= multi_X[, 7]),
      x8  <- as.numeric(0 <= multi_X[, 8]),
      x9  <- as.numeric(0 <= multi_X[, 9]),
      x10  <- as.numeric(0 <= multi_X[, 10])
    )
    
    # generate exposure
    logit_zp <- scenario$alpha0 + t(alpha %*% t(X))
    zp       <- exp(logit_zp) / (1 + exp(logit_zp))
    z        <- as.numeric(runif(scenario$n) < zp)
    
    # generate outcome
    logit_yp <- scenario$beta0 + betaz * z + t(beta %*% t(X))
    yp       <- exp(logit_yp) / (1 + exp(logit_yp))
    y        <- as.numeric(runif(scenario$n) < yp)
    
    # dataframe
    cov  <- data.frame(z, X)
    data <- data.frame(y, z, X)
    cov_1  <- cbind(rep(1, scenario$n), rep(1, scenario$n), X)
    cov_0  <- cbind(rep(1, scenario$n), rep(0, scenario$n), X)
    
    
    # formula 1
    formula_y     <- y ~ z + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 
    formula_z     <- z ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
    
    #----Counting 0 cell
    cell_count <- function(out, cov){
      sepa <- 0
      out <- factor(out, levels = c(0, 1))
      cov <- factor(cov, levels = c(0, 1))
      tab <- table(out, cov)
      if (any(tab == 0)) return(1) else return(0)
    }
    
    # ML
    detect_separation <- glm(formula = formula_y, data = data, family = binomial("logit"),
                             method = "detect_separation")
    for (j in 1 : length(detect_separation$coefficients)) {
      if(is.infinite(detect_separation$coefficients[j]) == TRUE){
        res_rsk_1$sepa[i] <- 1
        res_rsk_0$sepa[i] <- 1
        res_rsk_diff$sepa[i] <- 1
        res_rsk_log_ratio$sepa[i] <- 1
        break
      }
    }
    if(res_rsk_1$sepa[i] != 1){
      ML                      <- summary(glm(formula = formula_y, data = data, family = binomial(logit)))
      est_ML                  <- ML$coefficients[, "Estimate"]
      res_ML                  <- mrg_rsk_cal(est_ML, cov_1, cov_0)
      res_rsk_1$ML[i]         <- res_ML$rsk_1 
      res_rsk_0$ML[i]         <- res_ML$rsk_0 
      res_rsk_diff$ML[i]      <- res_ML$rsk_diff
      res_rsk_log_ratio$ML[i] <- res_ML$rsk_log_ratio
    }else{
      est_ML                  <- NA
      res_rsk_1$ML[i]         <- NA
      res_rsk_0$ML[i]         <- NA 
      res_rsk_diff$ML[i]      <- NA
      res_rsk_log_ratio$ML[i] <- NA
    }
    coverage_count_ML <- coverage_count(est = est_ML, cov_1 = cov_1, cov_0 = cov_0, covariance = ML$cov.unscaled,
                                        risk_1 = res_rsk_1$ML[i], risk_0 = res_rsk_0$ML[i], 
                                        risk_diff = res_rsk_diff$ML[i], risk_log_ratio = res_rsk_log_ratio$ML[i],
                                        n = scenario$n)
    coverage_count_rsk_1$ML[i]         <- coverage_count_ML$count_1
    coverage_count_rsk_0$ML[i]         <- coverage_count_ML$count_0
    coverage_count_rsk_diff$ML[i]      <- coverage_count_ML$count_diff
    coverage_count_rsk_log_ratio$ML[i] <- coverage_count_ML$count_log_ratio
    mese_count_rsk_1$ML[i]             <- coverage_count_ML$se_1
    mese_count_rsk_0$ML[i]             <- coverage_count_ML$se_0
    mese_count_rsk_diff$ML[i]          <- coverage_count_ML$se_diff
    mese_count_rsk_log_ratio$ML[i]     <- coverage_count_ML$se_log_ratio
    
    # no adj ML
    if(cell_count(y, z) == 1){
      est_no_adj_ML                  <- NA
      res_rsk_1$no_adj_ML[i]         <- NA
      res_rsk_0$no_adj_ML[i]         <- NA
      res_rsk_diff$no_adj_ML[i]      <- NA
      res_rsk_log_ratio$no_adj_ML[i] <- NA
    }else{
      no_adj_ML                      <- summary(glm(formula = y ~ z, data = data, family = binomial(logit)))
      est_no_adj_ML                  <- no_adj_ML$coefficients[, "Estimate"]
      res_no_adj_ML                  <- mrg_rsk_cal(est_no_adj_ML, cov_1, cov_0)
      res_rsk_1$no_adj_ML[i]         <- res_no_adj_ML$rsk_1
      res_rsk_0$no_adj_ML[i]         <- res_no_adj_ML$rsk_0
      res_rsk_diff$no_adj_ML[i]      <- res_no_adj_ML$rsk_diff
      res_rsk_log_ratio$no_adj_ML[i] <- res_no_adj_ML$rsk_log_ratio
    }    
    coverage_count_no_adj_ML <- coverage_count(est = est_no_adj_ML, cov_1 = cov_1[, 1 : 2], cov_0 = cov_0[, 1 : 2], covariance = no_adj_ML$cov.unscaled,
                                               risk_1 = res_rsk_1$no_adj_ML[i], risk_0 = res_rsk_0$no_adj_ML[i], 
                                               risk_diff = res_rsk_diff$no_adj_ML[i], risk_log_ratio = res_rsk_log_ratio$no_adj_ML[i],
                                               n = scenario$n)
    coverage_count_rsk_1$no_adj_ML[i]         <- coverage_count_no_adj_ML$count_1
    coverage_count_rsk_0$no_adj_ML[i]         <- coverage_count_no_adj_ML$count_0
    coverage_count_rsk_diff$no_adj_ML[i]      <- coverage_count_no_adj_ML$count_diff
    coverage_count_rsk_log_ratio$no_adj_ML[i] <- coverage_count_no_adj_ML$count_log_ratio
    mese_count_rsk_1$no_adj_ML[i]             <- coverage_count_no_adj_ML$se_1
    mese_count_rsk_0$no_adj_ML[i]             <- coverage_count_no_adj_ML$se_0
    mese_count_rsk_diff$no_adj_ML[i]          <- coverage_count_no_adj_ML$se_diff
    mese_count_rsk_log_ratio$no_adj_ML[i]     <- coverage_count_no_adj_ML$se_log_ratio
    
    
    # FML
    FML                      <- logistf(formula = formula_y, data = data, firth = T, pl = F)
    est_FML                  <- FML$coefficients
    res_FML                  <- mrg_rsk_cal(est_FML, cov_1, cov_0)
    res_rsk_1$FML[i]         <- res_FML$rsk_1 
    res_rsk_0$FML[i]         <- res_FML$rsk_0 
    res_rsk_diff$FML[i]      <- res_FML$rsk_diff
    res_rsk_log_ratio$FML[i] <- res_FML$rsk_log_ratio
    coverage_count_FML <- coverage_count(est = est_FML, cov_1 = cov_1, cov_0 = cov_0, covariance = FML$var,
                                         risk_1 = res_rsk_1$FML[i], risk_0 = res_rsk_0$FML[i], 
                                         risk_diff = res_rsk_diff$FML[i], risk_log_ratio = res_rsk_log_ratio$FML[i], 
                                         n = scenario$n)
    coverage_count_rsk_1$FML[i]         <- coverage_count_FML$count_1
    coverage_count_rsk_0$FML[i]         <- coverage_count_FML$count_0
    coverage_count_rsk_diff$FML[i]      <- coverage_count_FML$count_diff
    coverage_count_rsk_log_ratio$FML[i] <- coverage_count_FML$count_log_ratio
    mese_count_rsk_1$FML[i]             <- coverage_count_FML$se_1
    mese_count_rsk_0$FML[i]             <- coverage_count_FML$se_0
    mese_count_rsk_diff$FML[i]          <- coverage_count_FML$se_diff
    mese_count_rsk_log_ratio$FML[i]     <- coverage_count_FML$se_log_ratio
    
    # FLIC
    FLIC.fit1 <- FML
    FLIC.lp <- FLIC.fit1$linear.predictors - FLIC.fit1$coef[1]
    FLIC.fit2 <- glm(y ~ 1, family = binomial(link=logit), offset = FLIC.lp)
    W <- diag(FLIC.fit2$fitted.values * (1 - FLIC.fit2$fitted.values))
    FLIC.var <- solve(t(cbind(1,z,X)) %*% W %*% cbind(1,z,X))
    ic <- FLIC.fit2$coef
    FLIC <- list()
    FLIC$coefficients <- c(ic, FLIC.fit1$coef[-1])
    est_FLIC <- c()
    est_FLIC <- FLIC$coefficients
    res_FLIC                  <- mrg_rsk_cal(est_FLIC, cov_1, cov_0)
    res_rsk_1$FLIC[i]         <- res_FLIC$rsk_1 
    res_rsk_0$FLIC[i]         <- res_FLIC$rsk_0 
    res_rsk_diff$FLIC[i]      <- res_FLIC$rsk_diff
    res_rsk_log_ratio$FLIC[i] <- res_FLIC$rsk_log_ratio
    coverage_count_FLIC <- coverage_count(est = est_FLIC, cov_1 = cov_1, cov_0 = cov_0, covariance = FLIC.var,
                                          risk_1 = res_rsk_1$FLIC[i], risk_0 = res_rsk_0$FLIC[i], 
                                          risk_diff = res_rsk_diff$FLIC[i], risk_log_ratio = res_rsk_log_ratio$FLIC[i],
                                          n = scenario$n)
    coverage_count_rsk_1$FLIC[i]         <- coverage_count_FLIC$count_1
    coverage_count_rsk_0$FLIC[i]         <- coverage_count_FLIC$count_0
    coverage_count_rsk_diff$FLIC[i]      <- coverage_count_FLIC$count_diff
    coverage_count_rsk_log_ratio$FLIC[i] <- coverage_count_FLIC$count_log_ratio
    mese_count_rsk_1$FLIC[i]             <- coverage_count_FLIC$se_1
    mese_count_rsk_0$FLIC[i]             <- coverage_count_FLIC$se_0
    mese_count_rsk_diff$FLIC[i]          <- coverage_count_FLIC$se_diff
    mese_count_rsk_log_ratio$FLIC[i]     <- coverage_count_FLIC$se_log_ratio
    
    # FLAC
    FLAC.fit1 <- FML
    FLAC.pseudo <- c(rep(0, length(y)), rep(1,  2 * length(y)))
    FLAC.neww <- c(rep(1, length(y)), FLAC.fit1$hat/2, FLAC.fit1$hat/2)
    FLAC.fit2 <- logistf (formula = c(y, y, 1 - y) ~ c(z, z, z) + c(x1, x1, x1) + c(x2, x2, x2) + c(x3, x3, x3) + c(x4, x4, x4) + c(x5, x5, x5) +
                            c(x6, x6, x6) + c(x7, x7, x7) + c(x8, x8, x8) + c(x9, x9, x9) + c(x10, x10, x10) + FLAC.pseudo, 
                          weights = FLAC.neww, family = binomial(logit), firth = FALSE, pl = F)
    FLAC <- list()
    FLAC$coefficients <- FLAC.fit2$coefficients[which("FLAC.pseudo"!= names(FLAC.fit2$coefficients) )]
    est_FLAC <- c()
    est_FLAC <- FLAC$coefficients
    res_FLAC                  <- mrg_rsk_cal(est_FLAC, cov_1, cov_0)
    res_rsk_1$FLAC[i]         <- res_FLAC$rsk_1 
    res_rsk_0$FLAC[i]         <- res_FLAC$rsk_0 
    res_rsk_diff$FLAC[i]      <- res_FLAC$rsk_diff
    res_rsk_log_ratio$FLAC[i] <- res_FLAC$rsk_log_ratio
    coverage_count_FLAC <- coverage_count(est = est_FLAC, cov_1 = cov_1, cov_0 = cov_0, covariance = FLAC.fit2$var[1 : length(est_FLAC), 1 : length(est_FLAC)],
                                          risk_1 = res_rsk_1$FLAC[i], risk_0 = res_rsk_0$FLAC[i], 
                                          risk_diff = res_rsk_diff$FLAC[i], risk_log_ratio = res_rsk_log_ratio$FLAC[i],
                                          n = scenario$n)
    coverage_count_rsk_1$FLAC[i]         <- coverage_count_FLAC$count_1
    coverage_count_rsk_0$FLAC[i]         <- coverage_count_FLAC$count_0
    coverage_count_rsk_diff$FLAC[i]      <- coverage_count_FLAC$count_diff
    coverage_count_rsk_log_ratio$FLAC[i] <- coverage_count_FLAC$count_log_ratio
    mese_count_rsk_1$FLAC[i]             <- coverage_count_FLAC$se_1
    mese_count_rsk_0$FLAC[i]             <- coverage_count_FLAC$se_0
    mese_count_rsk_diff$FLAC[i]          <- coverage_count_FLAC$se_diff
    mese_count_rsk_log_ratio$FLAC[i]     <- coverage_count_FLAC$se_log_ratio
    
    # ps_adj & IPW
    if(cell_count(y, z) == 1){
      res_rsk_1$ps_adj[i]         <- NA 
      res_rsk_0$ps_adj[i]         <- NA 
      res_rsk_diff$ps_adj[i]      <- NA
      res_rsk_log_ratio$ps_adj[i] <- NA
      res_rsk_1$IPW[i]         <- NA
      res_rsk_0$IPW[i]         <- NA
      res_rsk_diff$IPW[i]      <- NA
      res_rsk_log_ratio$IPW[i] <- NA
    }else{
      ps <- fitted.values(glm(formula = formula_z, data = data, family = binomial(logit)))
      cov_ps_1 <- cbind(rep(1, scenario$n), rep(1, scenario$n), ps)
      cov_ps_0 <- cbind(rep(1, scenario$n), rep(0, scenario$n), ps)
      ps_adj <- summary(glm(formula = y ~ z + ps, data = data, family = binomial(logit)))
      est_ps_adj                  <- ps_adj$coefficients[, "Estimate"]
      res_ps_adj                  <- mrg_rsk_cal(est_ps_adj, cov_ps_1, cov_ps_0)
      res_rsk_1$ps_adj[i]         <- res_ps_adj$rsk_1 
      res_rsk_0$ps_adj[i]         <- res_ps_adj$rsk_0 
      res_rsk_diff$ps_adj[i]      <- res_ps_adj$rsk_diff
      res_rsk_log_ratio$ps_adj[i] <- res_ps_adj$rsk_log_ratio
      res_rsk_1$IPW[i]         <- (sum(data$z * y / ps)) / (sum(data$z / ps))
      res_rsk_0$IPW[i]         <- (sum((1 - data$z) * y / (1 - ps))) / (sum((1 - data$z) / (1 - ps)))
      res_rsk_diff$IPW[i]      <- res_rsk_1$IPW[i] - res_rsk_0$IPW[i]
      res_rsk_log_ratio$IPW[i] <- log(res_rsk_1$IPW[i] / res_rsk_0$IPW[i])
      IPW <- rep(NA, scenario$n)
      for (j in 1 : scenario$n) {
        IPW[j] <-((z[j] / ps[j]) + ((1 - z[j])/(1 - ps[j])))
      }
      data$id <- 1 : scenario$n
      res <- summary(geeglm(formula = y ~ z, id = id, data = data, family = binomial,
                            weights = IPW, std.err = "san.se"))
      est_IPW <- res$coefficients[, "Estimate"]
    }
    coverage_count_ps_adj <- coverage_count(est = est_ps_adj, cov_1 = cov_ps_1, cov_0 = cov_ps_0, covariance = ps_adj$cov.unscaled,
                                            risk_1 = res_rsk_1$ps_adj[i], risk_0 = res_rsk_0$ps_adj[i], 
                                            risk_diff = res_rsk_diff$ps_adj[i], risk_log_ratio = res_rsk_log_ratio$ps_adj[i],
                                            n = scenario$n)
    coverage_count_rsk_1$ps_adj[i]         <- coverage_count_ps_adj$count_1
    coverage_count_rsk_0$ps_adj[i]         <- coverage_count_ps_adj$count_0
    coverage_count_rsk_diff$ps_adj[i]      <- coverage_count_ps_adj$count_diff
    coverage_count_rsk_log_ratio$ps_adj[i] <- coverage_count_ps_adj$count_log_ratio
    mese_count_rsk_1$ps_adj[i]             <- coverage_count_ps_adj$se_1
    mese_count_rsk_0$ps_adj[i]             <- coverage_count_ps_adj$se_0
    mese_count_rsk_diff$ps_adj[i]          <- coverage_count_ps_adj$se_diff
    mese_count_rsk_log_ratio$ps_adj[i]     <- coverage_count_ps_adj$se_log_ratio
    coverage_count_IPW <- coverage_count(est = est_IPW, cov_1 = cov_1[, 1 : 2], cov_0 = cov_0[, 1 : 2], covariance = res$cov.unscaled,
                                         risk_1 = res_rsk_1$IPW[i], risk_0 = res_rsk_0$IPW[i], 
                                         risk_diff = res_rsk_diff$IPW[i], risk_log_ratio = res_rsk_log_ratio$IPW[i],
                                         n = scenario$n)
    coverage_count_rsk_1$IPW[i]         <- coverage_count_IPW$count_1
    coverage_count_rsk_0$IPW[i]         <- coverage_count_IPW$count_0
    coverage_count_rsk_diff$IPW[i]      <- coverage_count_IPW$count_diff
    coverage_count_rsk_log_ratio$IPW[i] <- coverage_count_IPW$count_log_ratio
    mese_count_rsk_1$IPW[i]             <- coverage_count_IPW$se_1
    mese_count_rsk_0$IPW[i]             <- coverage_count_IPW$se_0
    mese_count_rsk_diff$IPW[i]          <- coverage_count_IPW$se_diff
    mese_count_rsk_log_ratio$IPW[i]     <- coverage_count_IPW$se_log_ratio
    
    # # progress check
    # cat("\r", i)
  }
  
  #----Evaluation
  eval <- function(rsk_1, rsk_0, rsk_diff, rsk_log_ratio){
    na_rsk_1           <- 0
    na_rsk_0           <- 0
    na_rsk_diff        <- 0
    na_rsk_log_ratio   <- 0
    for (i in 1 : scenario$N_sim) {
      if(is.na(rsk_1[i]) == T || is.na(rsk_0[i]) == T){
        na_rsk_1 <- na_rsk_1 + 1
        na_rsk_0 <- na_rsk_0 + 1
        na_rsk_diff <- na_rsk_diff + 1
        na_rsk_log_ratio <- na_rsk_log_ratio + 1
        rsk_1[i] <- NA
        rsk_0[i] <- NA
        rsk_diff[i] <- NA
        rsk_log_ratio[i] <- NA
      }
    }
    mean_rsk_1         <- mean(rsk_1, na.rm = T)
    mean_rsk_0         <- mean(rsk_0, na.rm = T)
    mean_rsk_diff      <- mean(rsk_diff, na.rm = T)
    mean_rsk_log_ratio <- mean(rsk_log_ratio, na.rm = T)
    bias_rsk_1         <- mean(rsk_1, na.rm = T) - mu1
    bias_rsk_0         <- mean(rsk_0, na.rm = T) - mu0
    bias_rsk_diff      <- mean(rsk_diff, na.rm = T) - diff
    bias_rsk_log_ratio <- mean(rsk_log_ratio, na.rm = T) - log_ratio
    MCSE_rsk_1         <- sqrt(sum((rsk_1 - mean_rsk_1) ^ 2, na.rm = T) / (scenario$N_sim - 1 - na_rsk_1))
    MCSE_rsk_0         <- sqrt(sum((rsk_0 - mean_rsk_0) ^ 2, na.rm = T) / (scenario$N_sim - 1 - na_rsk_0))
    MCSE_rsk_diff      <- sqrt(sum((rsk_diff - mean_rsk_diff) ^ 2, na.rm = T) / (scenario$N_sim - 1 - na_rsk_diff))
    MCSE_rsk_log_ratio <- sqrt(sum((rsk_log_ratio - mean_rsk_log_ratio) ^ 2, na.rm = T) / (scenario$N_sim - 1 - na_rsk_log_ratio))
    rel_bias_rsk_1     <- bias_rsk_1 / mu1
    rel_bias_rsk_0     <- bias_rsk_0 / mu0
    rel_bias_rsk_diff  <- bias_rsk_diff / diff
    rel_bias_rsk_log_ratio <- bias_rsk_log_ratio / log_ratio 
    df <- data.frame(na_rsk_1 = na_rsk_1, na_rsk_0 = na_rsk_0, na_rsk_diff = na_rsk_diff, na_rsk_log_ratio = na_rsk_log_ratio, 
                     bias_rsk_1 = bias_rsk_1, bias_rsk_0 = bias_rsk_0, bias_rsk_diff = bias_rsk_diff, bias_rsk_log_ratio = bias_rsk_log_ratio, 
                     mean_rsk_1 = mean_rsk_1, mean_rsk_0 = mean_rsk_0, mean_rsk_diff = mean_rsk_diff, mean_rsk_log_ratio = mean_rsk_log_ratio, 
                     MCSE_rsk_1 = MCSE_rsk_1, MCSE_rsk_0 = MCSE_rsk_0, MCSE_rsk_diff = MCSE_rsk_diff, MCSE_rsk_log_ratio = MCSE_rsk_log_ratio,
                     rel_bias_rsk_1 = rel_bias_rsk_1, rel_bias_rsk_0 = rel_bias_rsk_0, rel_bias_rsk_diff = rel_bias_rsk_diff, rel_bias_rsk_log_ratio = rel_bias_rsk_log_ratio)
    return(df)
  }
  eval_ML         <- eval(res_rsk_1$ML, res_rsk_0$ML, res_rsk_diff$ML, res_rsk_log_ratio$ML)
  eval_no_adj_ML  <- eval(res_rsk_1$no_adj_ML, res_rsk_0$no_adj_ML, res_rsk_diff$no_adj_ML, res_rsk_log_ratio$no_adj_ML)
  eval_FML        <- eval(res_rsk_1$FML, res_rsk_0$FML, res_rsk_diff$FML, res_rsk_log_ratio$FML)
  eval_FLIC       <- eval(res_rsk_1$FLIC, res_rsk_0$FLIC, res_rsk_diff$FLIC, res_rsk_log_ratio$FLIC)
  eval_FLAC       <- eval(res_rsk_1$FLAC, res_rsk_0$FLAC, res_rsk_diff$FLAC, res_rsk_log_ratio$FLAC)
  eval_ps_adj     <- eval(res_rsk_1$ps_adj, res_rsk_0$ps_adj, res_rsk_diff$ps_adj, res_rsk_log_ratio$ps_adj)
  eval_IPW        <- eval(res_rsk_1$IPW, res_rsk_0$IPW, res_rsk_diff$IPW, res_rsk_log_ratio$IPW)
  
  #----Bias
  bias_rsk_1 <- data.frame(n = scenario$n,
                           measure = "mu1",
                           true_value = mu1,
                           Indicator = "bias",
                           ML = eval_ML$bias_rsk_1,
                           no_adj_ML = eval_no_adj_ML$bias_rsk_1,
                           FML = eval_FML$bias_rsk_1, 
                           FLIC = eval_FLIC$bias_rsk_1,
                           FLAC = eval_FLAC$bias_rsk_1,
                           ps_adj = eval_ps_adj$bias_rsk_1,
                           IPW = eval_IPW$bias_rsk_1)
  bias_rsk_0 <- data.frame(n = scenario$n,
                           measure = "mu0",
                           true_value = mu0,
                           Indicator = "bias",
                           ML = eval_ML$bias_rsk_0,
                           no_adj_ML = eval_no_adj_ML$bias_rsk_0,
                           FML = eval_FML$bias_rsk_0, 
                           FLIC = eval_FLIC$bias_rsk_0,
                           FLAC = eval_FLAC$bias_rsk_0,
                           ps_adj = eval_ps_adj$bias_rsk_0,
                           IPW = eval_IPW$bias_rsk_0)
  bias_rsk_diff <- data.frame(n = scenario$n,
                              measure = "diff",
                              true_value = diff,
                              Indicator = "bias",
                              ML = eval_ML$bias_rsk_diff,
                              no_adj_ML = eval_no_adj_ML$bias_rsk_diff,
                              FML = eval_FML$bias_rsk_diff, 
                              FLIC = eval_FLIC$bias_rsk_diff,
                              FLAC = eval_FLAC$bias_rsk_diff,
                              ps_adj = eval_ps_adj$bias_rsk_diff,
                              IPW = eval_IPW$bias_rsk_diff)
  bias_rsk_log_ratio <- data.frame(n = scenario$n,
                                   measure = "log_ratio",
                                   true_value = log_ratio,
                                   Indicator = "bias",
                                   ML = eval_ML$bias_rsk_log_ratio,
                                   no_adj_ML = eval_no_adj_ML$bias_rsk_log_ratio,
                                   FML = eval_FML$bias_rsk_log_ratio, 
                                   FLIC = eval_FLIC$bias_rsk_log_ratio,
                                   FLAC = eval_FLAC$bias_rsk_log_ratio,
                                   ps_adj = eval_ps_adj$bias_rsk_log_ratio,
                                   IPW = eval_IPW$bias_rsk_log_ratio)
  
  #----MCSE
  MCSE_rsk_1 <- data.frame(n = scenario$n,
                           measure = "mu1",
                           true_value = mu1,
                           Indicator = "MCSE",
                           ML = eval_ML$MCSE_rsk_1,
                           no_adj_ML = eval_no_adj_ML$MCSE_rsk_1,
                           FML = eval_FML$MCSE_rsk_1, 
                           FLIC = eval_FLIC$MCSE_rsk_1,
                           FLAC = eval_FLAC$MCSE_rsk_1,
                           ps_adj = eval_ps_adj$MCSE_rsk_1,
                           IPW = eval_IPW$MCSE_rsk_1)
  MCSE_rsk_0 <- data.frame(n = scenario$n,
                           measure = "mu0",
                           true_value = mu0,
                           Indicator = "MCSE",
                           ML = eval_ML$MCSE_rsk_0,
                           no_adj_ML = eval_no_adj_ML$MCSE_rsk_0,
                           FML = eval_FML$MCSE_rsk_0, 
                           FLIC = eval_FLIC$MCSE_rsk_0,
                           FLAC = eval_FLAC$MCSE_rsk_0,
                           ps_adj = eval_ps_adj$MCSE_rsk_0,
                           IPW = eval_IPW$MCSE_rsk_0)
  MCSE_rsk_diff <- data.frame(n = scenario$n,
                              measure = "diff",
                              true_value = diff,
                              Indicator = "MCSE",
                              ML = eval_ML$MCSE_rsk_diff,
                              no_adj_ML = eval_no_adj_ML$MCSE_rsk_diff,
                              FML = eval_FML$MCSE_rsk_diff, 
                              FLIC = eval_FLIC$MCSE_rsk_diff,
                              FLAC = eval_FLAC$MCSE_rsk_diff,
                              ps_adj = eval_ps_adj$MCSE_rsk_diff,
                              IPW = eval_IPW$MCSE_rsk_diff)
  MCSE_rsk_log_ratio <- data.frame(n = scenario$n,
                                   measure = "log_ratio",
                                   true_value = log_ratio,
                                   Indicator = "MCSE",
                                   ML = eval_ML$MCSE_rsk_log_ratio,
                                   no_adj_ML = eval_no_adj_ML$MCSE_rsk_log_ratio,
                                   FML = eval_FML$MCSE_rsk_log_ratio, 
                                   FLIC = eval_FLIC$MCSE_rsk_log_ratio,
                                   FLAC = eval_FLAC$MCSE_rsk_log_ratio,
                                   ps_adj = eval_ps_adj$MCSE_rsk_log_ratio,
                                   IPW = eval_IPW$MCSE_rsk_log_ratio)
  #----relative_bias
  rel_bias_rsk_1 <- data.frame(n = scenario$n,
                               measure = "mu1",
                               true_value = mu1,
                               Indicator = "relative_bias",
                               ML = eval_ML$rel_bias_rsk_1,
                               no_adj_ML = eval_no_adj_ML$rel_bias_rsk_1,
                               FML = eval_FML$rel_bias_rsk_1, 
                               FLIC = eval_FLIC$rel_bias_rsk_1,
                               FLAC = eval_FLAC$rel_bias_rsk_1,
                               ps_adj = eval_ps_adj$rel_bias_rsk_1,
                               IPW = eval_IPW$rel_bias_rsk_1)
  rel_bias_rsk_0 <- data.frame(n = scenario$n,
                               measure = "mu0",
                               true_value = mu0,
                               Indicator = "relative_bias",
                               ML = eval_ML$rel_bias_rsk_0,
                               no_adj_ML = eval_no_adj_ML$rel_bias_rsk_0,
                               FML = eval_FML$rel_bias_rsk_0, 
                               FLIC = eval_FLIC$rel_bias_rsk_0,
                               FLAC = eval_FLAC$rel_bias_rsk_0,
                               ps_adj = eval_ps_adj$rel_bias_rsk_0,
                               IPW = eval_IPW$rel_bias_rsk_0)
  rel_bias_rsk_diff <- data.frame(n = scenario$n,
                                  measure = "diff",
                                  true_value = diff,
                                  Indicator = "relative_bias",
                                  ML = eval_ML$rel_bias_rsk_diff,
                                  no_adj_ML = eval_no_adj_ML$rel_bias_rsk_diff,
                                  FML = eval_FML$rel_bias_rsk_diff, 
                                  FLIC = eval_FLIC$rel_bias_rsk_diff,
                                  FLAC = eval_FLAC$rel_bias_rsk_diff,
                                  ps_adj = eval_ps_adj$rel_bias_rsk_diff,
                                  IPW = eval_IPW$rel_bias_rsk_diff)
  rel_bias_rsk_log_ratio <- data.frame(n = scenario$n,
                                       measure = "log_ratio",
                                       true_value = log_ratio,
                                       Indicator = "relative_bias",
                                       ML = eval_ML$rel_bias_rsk_log_ratio,
                                       no_adj_ML = eval_no_adj_ML$rel_bias_rsk_log_ratio,
                                       FML = eval_FML$rel_bias_rsk_log_ratio, 
                                       FLIC = eval_FLIC$rel_bias_rsk_log_ratio,
                                       FLAC = eval_FLAC$rel_bias_rsk_log_ratio,
                                       ps_adj = eval_ps_adj$rel_bias_rsk_log_ratio,
                                       IPW = eval_IPW$rel_bias_rsk_log_ratio)
  #----mese
  mese_rsk_1 <- data.frame(n = scenario$n,
                           measure = "mu1",
                           Indicator = "mese",
                           ML = mean(mese_count_rsk_1$ML, na.rm = T),
                           no_adj_ML = mean(mese_count_rsk_1$no_adj_ML, na.rm = T),
                           FML = mean(mese_count_rsk_1$FML, na.rm = T), 
                           FLIC = mean(mese_count_rsk_1$FLIC, na.rm = T),
                           FLAC = mean(mese_count_rsk_1$FLAC, na.rm = T),
                           ps_adj = mean(mese_count_rsk_1$ps_adj, na.rm = T),
                           IPW = mean(mese_count_rsk_1$IPW, na.rm = T))
  mese_rsk_0 <- data.frame(n = scenario$n,
                           measure = "mu0",
                           Indicator = "mese",
                           ML = mean(mese_count_rsk_0$ML, na.rm = T),
                           no_adj_ML = mean(mese_count_rsk_0$no_adj_ML, na.rm = T),
                           FML = mean(mese_count_rsk_0$FML, na.rm = T), 
                           FLIC = mean(mese_count_rsk_0$FLIC, na.rm = T),
                           FLAC = mean(mese_count_rsk_0$FLAC, na.rm = T),
                           ps_adj = mean(mese_count_rsk_0$ps_adj, na.rm = T),
                           IPW = mean(mese_count_rsk_0$IPW, na.rm = T))
  mese_rsk_diff <- data.frame(n = scenario$n,
                              measure = "diff",
                              Indicator = "mese",
                              ML = mean(mese_count_rsk_diff$ML, na.rm = T),
                              no_adj_ML = mean(mese_count_rsk_diff$no_adj_ML, na.rm = T),
                              FML = mean(mese_count_rsk_diff$FML, na.rm = T), 
                              FLIC = mean(mese_count_rsk_diff$FLIC, na.rm = T),
                              FLAC = mean(mese_count_rsk_diff$FLAC, na.rm = T),
                              ps_adj = mean(mese_count_rsk_diff$ps_adj, na.rm = T),
                              IPW = mean(mese_count_rsk_diff$IPW, na.rm = T))
  mese_rsk_log_ratio <- data.frame(n = scenario$n,
                                   measure = "log_ratio",
                                   Indicator = "mese",
                                   ML = mean(mese_count_rsk_log_ratio$ML, na.rm = T),
                                   no_adj_ML = mean(mese_count_rsk_log_ratio$no_adj_ML, na.rm = T),
                                   FML = mean(mese_count_rsk_log_ratio$FML, na.rm = T), 
                                   FLIC = mean(mese_count_rsk_log_ratio$FLIC, na.rm = T),
                                   FLAC = mean(mese_count_rsk_log_ratio$FLAC, na.rm = T),
                                   ps_adj = mean(mese_count_rsk_log_ratio$ps_adj, na.rm = T),
                                   IPW = mean(mese_count_rsk_log_ratio$IPW, na.rm = T))
  #----coverage
  coverage_rsk_1 <- data.frame(n = scenario$n,
                               measure = "mu1",
                               Indicator = "coverage",
                               ML = mean(coverage_count_rsk_1$ML),
                               no_adj_ML = mean(coverage_count_rsk_1$no_adj_ML),
                               FML = mean(coverage_count_rsk_1$FML), 
                               FLIC = mean(coverage_count_rsk_1$FLIC),
                               FLAC = mean(coverage_count_rsk_1$FLAC),
                               ps_adj = mean(coverage_count_rsk_1$ps_adj),
                               IPW = mean(coverage_count_rsk_1$IPW))
  coverage_rsk_0 <- data.frame(n = scenario$n,
                               measure = "mu0",
                               Indicator = "coverage",
                               ML = mean(coverage_count_rsk_0$ML),
                               no_adj_ML = mean(coverage_count_rsk_0$no_adj_ML),
                               FML = mean(coverage_count_rsk_0$FML), 
                               FLIC = mean(coverage_count_rsk_0$FLIC),
                               FLAC = mean(coverage_count_rsk_0$FLAC),
                               ps_adj = mean(coverage_count_rsk_0$ps_adj),
                               IPW = mean(coverage_count_rsk_0$IPW))
  coverage_rsk_diff <- data.frame(n = scenario$n,
                                  measure = "diff",
                                  Indicator = "coverage",
                                  ML = mean(coverage_count_rsk_diff$ML),
                                  no_adj_ML = mean(coverage_count_rsk_diff$no_adj_ML),
                                  FML = mean(coverage_count_rsk_diff$FML), 
                                  FLIC = mean(coverage_count_rsk_diff$FLIC),
                                  FLAC = mean(coverage_count_rsk_diff$FLAC),
                                  ps_adj = mean(coverage_count_rsk_diff$ps_adj),
                                  IPW = mean(coverage_count_rsk_diff$IPW))
  coverage_rsk_log_ratio <- data.frame(n = scenario$n,
                                       measure = "log_ratio",
                                       Indicator = "coverage",
                                       ML = mean(coverage_count_rsk_log_ratio$ML),
                                       no_adj_ML = mean(coverage_count_rsk_log_ratio$no_adj_ML),
                                       FML = mean(coverage_count_rsk_log_ratio$FML), 
                                       FLIC = mean(coverage_count_rsk_log_ratio$FLIC),
                                       FLAC = mean(coverage_count_rsk_log_ratio$FLAC),
                                       ps_adj = mean(coverage_count_rsk_log_ratio$ps_adj),
                                       IPW = mean(coverage_count_rsk_log_ratio$IPW))
  
  #----Ratio of between MCSE and mese
  ratio_SE_rsk_1 <- data.frame(n = scenario$n,
                               measure = "mu(1)",
                               ML = MCSE_rsk_1$ML / mese_rsk_1$ML,
                               no_adj_ML = MCSE_rsk_1$no_adj_ML / mese_rsk_1$no_adj_ML,
                               FML = MCSE_rsk_1$FML / mese_rsk_1$FML, 
                               FLIC = MCSE_rsk_1$FLIC / mese_rsk_1$FLIC,
                               FLAC = MCSE_rsk_1$FLAC / mese_rsk_1$FLAC,
                               ps_adj = MCSE_rsk_1$ps_adj / mese_rsk_1$ps_adj,
                               IPW = MCSE_rsk_1$IPW / mese_rsk_1$IPW)
  ratio_SE_rsk_0 <- data.frame(n = scenario$n,
                               measure = "mu(0)",
                               ML = MCSE_rsk_0$ML / mese_rsk_0$ML,
                               no_adj_ML = MCSE_rsk_0$no_adj_ML / mese_rsk_0$no_adj_ML,
                               FML = MCSE_rsk_0$FML / mese_rsk_0$FML, 
                               FLIC = MCSE_rsk_0$FLIC / mese_rsk_0$FLIC,
                               FLAC = MCSE_rsk_0$FLAC / mese_rsk_0$FLAC,
                               ps_adj = MCSE_rsk_0$ps_adj / mese_rsk_0$ps_adj,
                               IPW = MCSE_rsk_0$IPW / mese_rsk_0$IPW)
  ratio_SE_rsk_diff <- data.frame(n = scenario$n,
                                  measure = "RD",
                                  ML = MCSE_rsk_diff$ML / mese_rsk_diff$ML,
                                  no_adj_ML = MCSE_rsk_diff$no_adj_ML / mese_rsk_diff$no_adj_ML,
                                  FML = MCSE_rsk_diff$FML / mese_rsk_diff$FML, 
                                  FLIC = MCSE_rsk_diff$FLIC / mese_rsk_diff$FLIC,
                                  FLAC = MCSE_rsk_diff$FLAC / mese_rsk_diff$FLAC,
                                  ps_adj = MCSE_rsk_diff$ps_adj / mese_rsk_diff$ps_adj,
                                  IPW = MCSE_rsk_diff$IPW / mese_rsk_diff$IPW)
  ratio_SE_rsk_log_ratio <- data.frame(n = scenario$n,
                                       measure = "log(RR)",
                                       ML = MCSE_rsk_log_ratio$ML / mese_rsk_log_ratio$ML,
                                       no_adj_ML = MCSE_rsk_log_ratio$no_adj_ML / mese_rsk_log_ratio$no_adj_ML,
                                       FML = MCSE_rsk_log_ratio$FML / mese_rsk_log_ratio$FML, 
                                       FLIC = MCSE_rsk_log_ratio$FLIC / mese_rsk_log_ratio$FLIC,
                                       FLAC = MCSE_rsk_log_ratio$FLAC / mese_rsk_log_ratio$FLAC,
                                       ps_adj = MCSE_rsk_log_ratio$ps_adj / mese_rsk_log_ratio$ps_adj,
                                       IPW = MCSE_rsk_log_ratio$IPW / mese_rsk_log_ratio$IPW)
  # ---- file name
  file_name_est  <- sprintf("est-scenario-%d.csv", row_index)
  file_name_bias <- sprintf("bias-scenario-%d.csv", row_index)
  file_name_MCSE <- sprintf("MCSE-scenario-%d.csv", row_index)
  file_name_rb   <- sprintf("relative-bias-scenario-%d.csv", row_index)
  file_name_mese <- sprintf("mese-scenario-%d.csv", row_index)
  file_name_cove <- sprintf("coverage-scenario-%d.csv", row_index)
  file_name_rSE  <- sprintf("ratioSE-scenario-%d.csv", row_index)
  
  write.csv(est,          file = file.path(out_dir, file_name_est),  row.names = FALSE)
  write.csv(bias,         file = file.path(out_dir, file_name_bias), row.names = FALSE)
  write.csv(MCSE,         file = file.path(out_dir, file_name_MCSE), row.names = FALSE)
  write.csv(relative_bias,file = file.path(out_dir, file_name_rb),   row.names = FALSE)
  write.csv(mese,         file = file.path(out_dir, file_name_mese), row.names = FALSE)
  write.csv(coverage,     file = file.path(out_dir, file_name_cove), row.names = FALSE)
  write.csv(ratio_SE,     file = file.path(out_dir, file_name_rSE),  row.names = FALSE)
  
  return(list(
    est = file_name_est,
    bias = file_name_bias,
    MCSE = file_name_MCSE,
    relative_bias = file_name_rb,
    mese = file_name_mese,
    coverage = file_name_cove,
    ratio_SE = file_name_rSE
  ))
}

# ---- parallel execution ----
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c("sim_summary", "run_simulation", "data_dir", "out_dir"))
clusterEvalQ(cl, {
  library(MASS)
  library(logistf)
  library(geepack)
  library(detectseparation)
})

results <- pblapply(1:nrow(sim_summary), run_simulation, cl = cl)

stopCluster(cl)

print(results)
