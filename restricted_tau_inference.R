res.tau.hat_func <- function(X, observed.time, delta, t.star, alpha = 0.05) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  
  ## point estimation
  tau.hat <- res.tau.hat.est(X, observed.time, delta, t.star)
  
  ## var
  var3 <- res.tau.hat.var3.fixed(X, observed.time, delta, tau.hat, t.star)
  var2 <- res.tau.hat.var2.fixed(X, observed.time, delta, t.star)
  var.est <- (var3 - var2) / (N0 + N1)
  
  ## confidence interval
  ci.l <- tau.hat + qnorm(p = alpha / 2) * sqrt(var.est)
  ci.r <- tau.hat + qnorm(p = 1 - alpha / 2) * sqrt(var.est)
  
  return(list(tau.hat = tau.hat,
              var.est = var.est,
              ci = c(ci.l, ci.r)))
}

res.tau.hat.est <- function(X, observed.time, delta, t.star) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## restricted region indicator
  restricted.indicator <- ifelse(obs.min.time.mat <= t.star, 1, 0)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## stat
  U <- ifelse(orderable.indicator * restricted.indicator, conc / G0.val / G1.val, 0)
  tau.hat <- sum(U) / N0 / N1
  
  return(tau.hat)
}

res.tau.hat.var2.fixed <- function(X, observed.time, delta, t.star) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## var1 estimation
  eta_func <- function(u) {
    trunc.indicator <- (obs.min.time.mat >= u) & (t.star >= obs.min.time.mat)
    
    eta.val <- ifelse(orderable.indicator & trunc.indicator,
                      conc / G0.val / G1.val, 0)
    
    return(sum(eta.val) / N0 / N1)
  }
  
  eta.val <- sapply(X = observed.time, FUN = eta_func)
  R0.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.0), 2, sum)
  R1.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.1), 2, sum)
  var2.0 <- N0 * sum(ifelse((1 - delta) * (1 - X), (eta.val / R0.val) ^ 2, 0))
  var2.1 <- N1 * sum(ifelse((1 - delta) * (X), (eta.val / R1.val) ^ 2, 0))
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var2 <- var2.0 / p0.hat + var2.1 / p1.hat
  
  return(var2)
}

res.tau.hat.var3.fixed <- function(X, observed.time, delta, res.tau.hat, t.star) {
  n <- length(X)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## var3 estimation
  var3.01 <- 0
  for(j in 1:length(observed.time.1)) {
    min.index <- sapply(observed.time.0, "<=", observed.time.1[j])
    obs.min.time <- ifelse(min.index, observed.time.0, observed.time.1[j])
    orderable.indicator <- ifelse(min.index, delta.0, delta.1[j])
    orderable.indicator <- ifelse(observed.time.0 == observed.time.1[j] & delta.0 == delta.1[j],
                                  FALSE,
                                  orderable.indicator)
    restricted.region.indicator <- (obs.min.time <= t.star)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index, 1, -1)
    U <- ifelse(orderable.indicator & restricted.region.indicator,
                conc / G0.val / G1.val,
                0)
    expt <- mean(sapply(U, "*", U))
    var3.01 <- var3.01 + expt
  }
  sigma01.square <- var3.01 / N1 - res.tau.hat ^ 2
  
  var3.10 <- 0
  for(i in 1:length(observed.time.0)) {
    min.index <- sapply(observed.time.1, "<=", observed.time.0[i])
    obs.min.time <- ifelse(min.index, observed.time.1, observed.time.0[i])
    orderable.indicator <- ifelse(min.index, delta.1, delta.0[i])
    orderable.indicator <- ifelse(observed.time.1 == observed.time.0[i] & delta.1 == delta.0[i],
                                  FALSE,
                                  orderable.indicator)
    restricted.region.indicator <- (obs.min.time <= t.star)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index, -1, 1)
    U <- ifelse(orderable.indicator & restricted.region.indicator,
                conc / G0.val / G1.val,
                0)
    expt <- mean(sapply(U, "*", U))
    var3.10 <- var3.10 + expt
  }
  sigma10.square <- var3.10 / N0 - res.tau.hat ^ 2
  
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var3 <- sigma01.square / p1.hat + sigma10.square / p0.hat
  
  return(var3)
}

tau.hat.var2.fixed <- function(X, observed.time, delta) {
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  observed.time.0 <- observed.time[X == 0]
  delta.0 <- delta[X == 0]
  observed.time.1 <- observed.time[X == 1]
  delta.1 <- delta[X == 1]
  
  ## estiamte the survival functions of censoring variable
  G0.fit <- survfit(Surv(observed.time.0, 1 - delta.0) ~ 1)
  G1.fit <- survfit(Surv(observed.time.1, 1 - delta.1) ~ 1)
  G0.fit_func <- stepfun(G0.fit$time, c(1, G0.fit$surv))
  G1.fit_func <- stepfun(G1.fit$time, c(1, G1.fit$surv))
  
  ## all combinations
  obs.time.mat.0 <- matrix(observed.time.0, nrow = N0, ncol = N1)
  obs.time.mat.1 <- matrix(observed.time.1, nrow = N0, ncol = N1, byrow = TRUE)
  delta.mat.0 <- matrix(as.logical(delta.0), nrow = N0, ncol = N1)
  delta.mat.1 <- matrix(as.logical(delta.1), nrow = N0, ncol = N1, byrow = TRUE)
  
  ## tilde{Y}
  min.index <- obs.time.mat.0 <= obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index, delta.mat.0, delta.mat.1)
  orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 == delta.mat.1,
                                FALSE,
                                orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index, 1, -1)
  
  ## var1 estimation
  eta_func <- function(u) {
    trunc.indicator <- (obs.min.time.mat >= u)
    
    eta.val <- ifelse(orderable.indicator & trunc.indicator,
                      conc / G0.val / G1.val, 0)
    
    return(sum(eta.val) / N0 / N1)
  }
  
  eta.val <- sapply(X = observed.time, FUN = eta_func)
  R0.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.0), 2, sum)
  R1.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.1), 2, sum)
  var2.0 <- N0 * sum(ifelse((1 - delta) * (1 - X), (eta.val / R0.val) ^ 2, 0))
  var2.1 <- N1 * sum(ifelse((1 - delta) * (X), (eta.val / R1.val) ^ 2, 0))
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var2 <- var2.0 / p0.hat + var2.1 / p1.hat
  
  return(var2)
}
