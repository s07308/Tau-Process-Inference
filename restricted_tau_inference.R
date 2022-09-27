res.tau.process_func <- function(X, observed.time, delta, t.star) {
  tau.process <- numeric(length = length(t.star))
  
  for(k in seq_along(t.star)) {
    tau.process[k] <- res.tau.hat.est(X, observed.time, delta, t.star[k])
  }
  
  km.fit <- survfit(Surv(observed.time, delta) ~ X)
  
  par(mfrow = c(1, 2))
  plot(km.fit, lty = c(2, 1), main = "Kaplan-Meier Curves")
  legend("topright", legend = c("control", "treatment"), lty = c(2, 1))
  plot(t.star, tau.process, type = "b", main = "Tau process")
  abline(a = 0, b = 0, lty = 3)
  
  return(tau.process)
}

res.tau.hat_func <- function(X, observed.time, delta, t.star, alpha = 0.05) {
  n <- length(X)
  N0 <- sum(X == 0)
  N1 <- sum(X == 1)
  
  ## point estimation
  tau.hat <- res.tau.hat.est(X, observed.time, delta, t.star)
  
  ## var (fixed)
  var3.fixed <- res.tau.hat.var3.fixed(X, observed.time, delta, tau.hat, t.star)
  var2.fixed <- res.tau.hat.var2.fixed(X, observed.time, delta, t.star)
  var.est.fixed <- (var3.fixed - var2.fixed) / (N0 + N1)
  
  ## var (random)
  var3.random <- res.tau.hat.var3.random(X, observed.time, delta, tau.hat, t.star)
  var2.random <- res.tau.hat.var2.random(X, observed.time, delta, t.star)
  var1.random <- res.tau.hat.var1.random(X, tau.hat)
  var.est.random <- (var3.random - var2.random - var1.random) / n
  
  ## confidence interval (fixed)
  ci.l.fixed <- tau.hat + qnorm(p = alpha / 2) * sqrt(var.est.fixed)
  ci.r.fixed <- tau.hat + qnorm(p = 1 - alpha / 2) * sqrt(var.est.fixed)
  
  ## confidence interval (random)
  ci.l.random <- tau.hat + qnorm(p = alpha / 2) * sqrt(var.est.random)
  ci.r.random <- tau.hat + qnorm(p = 1 - alpha / 2) * sqrt(var.est.random)
  
  ## p-values (random)
  var3.null <- res.tau.hat.var3.random(X, observed.time, delta, 0, t.star)
  var.null <- (var3.null - var2.random) / n
  
  z.est <- tau.hat / sqrt(var.est.random)
  z.0 <- tau.hat / sqrt(var.null)
  p.est <- 1 - pchisq(z.est ^ 2, df = 1)
  p.0 <- 1 - pchisq(z.0 ^ 2, df = 1)
  
  ## the proportion being compared
  max.fail.0 <- max(observed.time[X == 0 & delta == 1])
  max.fail.1 <- max(observed.time[X == 1 & delta == 1])
  fail.cut <- min(c(max.fail.0, max.fail.1))
  km.fit <- survfit(Surv(observed.time, delta) ~ X)
  prop.compared <- ifelse(prod(summary(km.fit, times = fail.cut)$surv) == 0,
                          1, 1 - prod(summary(km.fit, times = t.star)$surv))
  
  
  return(list(tau.hat = tau.hat,
              var.est = var.est.random,
              ci.random = c(ci.l.random, ci.r.random),
              ci.fixed = c(ci.l.fixed, ci.r.fixed),
              p.est = p.est,
              p.0 = p.0, 
              prop.compared = prop.compared))
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
  min.index.0 <- obs.time.mat.0 < obs.time.mat.1
  min.index.1 <- obs.time.mat.1 < obs.time.mat.0
  obs.min.time.mat <- ifelse(min.index.0, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index.0, delta.mat.0, 0)
  orderable.indicator <- ifelse(min.index.1, delta.mat.1, orderable.indicator)
  # orderable.indicator <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 != delta.mat.1,
  #                               1,
  #                               orderable.indicator)
  
  ## restricted region indicator
  restricted.indicator <- ifelse(obs.min.time.mat <= t.star, 1, 0)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index.0, 1, 0)
  conc <- ifelse(min.index.1, -1, conc)
  # conc <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 > delta.mat.1,
  #                1, conc)
  # conc <- ifelse(obs.time.mat.0 == obs.time.mat.1 & delta.mat.0 < delta.mat.1,
  #                -1, conc)
  
  
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
  min.index.0 <- obs.time.mat.0 < obs.time.mat.1
  min.index.1 <- obs.time.mat.1 < obs.time.mat.0
  obs.min.time.mat <- ifelse(min.index.0, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index.0, delta.mat.0, 0)
  orderable.indicator <- ifelse(min.index.1, delta.mat.1, orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index.0, 1, 0)
  conc <- ifelse(min.index.1, -1, conc)
  
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
    min.index.0 <- sapply(observed.time.0, "<", observed.time.1[j])
    min.index.1 <- sapply(observed.time.0, ">", observed.time.1[j])
    obs.min.time <- ifelse(min.index.0, observed.time.0, observed.time.1[j])
    orderable.indicator <- ifelse(min.index.0, delta.0, 0)
    orderable.indicator <- ifelse(min.index.1, delta.1[j], orderable.indicator)
    restricted.region.indicator <- (obs.min.time <= t.star)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index.0, 1, 0)
    conc <- ifelse(min.index.1, -1, conc)
    U <- ifelse(orderable.indicator & restricted.region.indicator,
                conc / G0.val / G1.val,
                0)
    expt <- mean(sapply(U, "*", U))
    var3.01 <- var3.01 + expt
  }
  sigma01.square <- var3.01 / N1 - res.tau.hat ^ 2
  
  var3.10 <- 0
  for(i in 1:length(observed.time.0)) {
    min.index.1 <- sapply(observed.time.1, "<", observed.time.0[i])
    min.index.0 <- sapply(observed.time.1, ">", observed.time.0[i])
    obs.min.time <- ifelse(min.index.1, observed.time.1, observed.time.0[i])
    orderable.indicator <- ifelse(min.index.1, delta.1, 0)
    orderable.indicator <- ifelse(min.index.0, delta.0[i], orderable.indicator)
    restricted.region.indicator <- (obs.min.time <= t.star)
    G0.val <- sapply(obs.min.time, G0.fit_func)
    G1.val <- sapply(obs.min.time, G1.fit_func)
    conc <- ifelse(min.index.0, 1, 0)
    conc <- ifelse(min.index.1, -1, conc)
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

res.tau.hat.var1.random <- function(X, res.tau.hat) {
  ## var1 estimation
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var1 <- res.tau.hat ^ 2 * (p1.hat - p0.hat) ^ 2 / (p0.hat * p1.hat)
  
  return(var1)
}

res.tau.hat.var2.random <- function(X, observed.time, delta, t.star) {
  n <- length(observed.time)
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
  min.index.0 <- obs.time.mat.0 < obs.time.mat.1
  min.index.1 <- obs.time.mat.0 > obs.time.mat.1
  obs.min.time.mat <- ifelse(min.index.0, obs.time.mat.0, obs.time.mat.1)
  
  ## orderable indicator
  orderable.indicator <- ifelse(min.index.0, delta.mat.0, 0)
  orderable.indicator <- ifelse(min.index.1, delta.mat.1, orderable.indicator)
  
  ## G0.hat and G1.hat at tilde{Y}
  G0.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G0.fit_func)
  G1.val <- apply(X = obs.min.time.mat, MARGIN = 2, FUN = G1.fit_func)
  
  ## concordance and discordance
  conc <- ifelse(min.index.0, 1, 0)
  conc <- ifelse(min.index.1, -1, conc)
  
  ## var2 estimation
  kappa_func <- function(u) {
    trunc.indicator <- (obs.min.time.mat >= u) & (t.star >= obs.min.time.mat)
    
    kappa.val <- ifelse(orderable.indicator & trunc.indicator,
                        conc / G0.val / G1.val, 0)
    
    return(sum(kappa.val) / choose(n, 2))
  }
  
  kappa.val <- sapply(X = observed.time, FUN = kappa_func)
  R0.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.0), 2, sum)
  R1.val <- apply(sapply(X = observed.time, FUN = "<=", observed.time.1), 2, sum)
  var2.0 <- n * sum(ifelse((1 - delta) * (1 - X), (kappa.val / R0.val) ^ 2, 0))
  var2.1 <- n * sum(ifelse((1 - delta) * (X), (kappa.val / R1.val) ^ 2, 0))
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  var2 <- (var2.0 + var2.1) / (4 * p0.hat ^ 2 * p1.hat ^ 2)
  
  return(var2)
}

res.tau.hat.var3.random <- function(X, observed.time, delta, res.tau.hat, t.star) {
  n <- length(observed.time)
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
  var3 <- 0
  for(k in 1:length(X)) {
    if(X[k] == 0) {
      min.index.1 <- sapply(observed.time.1, "<", observed.time[k])
      min.index.0 <- sapply(observed.time.1, ">", observed.time[k])
      obs.min.time <- ifelse(min.index.1, observed.time.1, observed.time[k])
      orderable.indicator <- ifelse(min.index.1, delta.1, 0)
      orderable.indicator <- ifelse(min.index.0, delta[k], orderable.indicator)
      restricted.region.indicator <- (obs.min.time <= t.star)
      G0.val <- sapply(obs.min.time, G0.fit_func)
      G1.val <- sapply(obs.min.time, G1.fit_func)
      conc <- ifelse(min.index.1, -1, 0)
      conc <- ifelse(min.index.0, 1, conc)
      U <- ifelse(orderable.indicator & restricted.region.indicator, conc / G0.val / G1.val, 0)
      expt <- sum(sapply(U, "*", U)) / n ^ 2
      # expt <- sum(U) / n
      var3 <- var3 + expt
      # var3 <- var3 + expt ^ 2
    }
    
    if(X[k] == 1) {
      min.index.0 <- sapply(observed.time.0, "<", observed.time[k])
      min.index.1 <- sapply(observed.time.0, ">", observed.time[k])
      obs.min.time <- ifelse(min.index.0, observed.time.0, observed.time[k])
      orderable.indicator <- ifelse(min.index.0, delta.0, 0)
      orderable.indicator <- ifelse(min.index.1, delta[k], orderable.indicator)
      restricted.region.indicator <- (obs.min.time <= t.star)
      G0.val <- sapply(obs.min.time, G0.fit_func)
      G1.val <- sapply(obs.min.time, G1.fit_func)
      conc <- ifelse(min.index.0, 1, 0)
      conc <- ifelse(min.index.1, -1, conc)
      U <- ifelse(orderable.indicator & restricted.region.indicator, conc / G0.val / G1.val, 0)
      expt <- sum(sapply(U, "*", U)) / n ^ 2
      # expt <- sum(U) / n
      var3 <- var3 + expt
      # var3 <- var3 + expt ^ 2
    }
  }
  
  p1.hat <- mean(X == 1)
  p0.hat <- 1 - p1.hat
  theta1.square <- var3 / n - (2 * p0.hat * p1.hat * res.tau.hat) ^ 2
  var3 <- theta1.square / (p0.hat * p1.hat) ^ 2
  
  return(var3)
}

imputed.tau.hat_func2 <- function(X, observed.time, delta, t.min, t.max, dist.0, dist.1) {
  df <- data.frame(time = observed.time, censor = delta, arm = X)
  df.0 <- df[df$X == 0, ]
  df.1 <- df[df$X == 1, ]
  param.fit.0 <- fit_data(data = df.0, dist = dist.0)
  param.fit.1 <- fit_data(data = df.1, dist = dist.1)
  d0 <- dens_fun(dist = dist.0, param.fit = param.fit.0)
  s0 <- surv_fun(dist = dist.0, param.fit = param.fit.0)
  d1 <- dens_fun(dist = dist.1, param.fit = param.fit.1)
  s1 <- surv_fun(dist = dist.1, param.fit = param.fit.1)
  
  imputed.tail <- integrate(f = function(t) {
    d0(t) * s1(t) - d1(t) * s0(t)
  },
  lower = t.min,
  upper = t.max) 
  
  res.tau.hat <- res.tau.hat.est(X = arm, observed.time = surv.time, delta = event, t.star = t.min)
  
  return(res.tau.hat + imputed.tail$value)  
}

dens_fun <- function(dist, param.fit) {
  if(dist == "weibull") {
    return(
      function(x) dweibull(x = x, shape = param.fit$estimate["shape"], scale = param.fit$estimate["scale"])
    )
  }
  
  if(dist == "exp") {
    return(
      function(x) dexp(x = x, rate = param.fit$estimate["rate"])
    )
  }
  
  if(dist == "lnorm") {
    return(
      function(x) dlnorm(x = x, meanlog = param.fit$estimate["meanlog"], sdlog = param.fit$estimate["sdlog"])
    )
  }
  
  if(dist == "logis") {
    return(
      function(x) dlogis(x = x, location = param.fit$estimate["location"], scale = param.fit$estimate["scale"])
    )
  }
}

surv_fun <- function(dist, param.fit) {
  if(dist == "weibull") {
    return(
      function(x) pweibull(q = x,
                           shape = param.fit$estimate["shape"],
                           scale = param.fit$estimate["scale"], 
                           lower.tail = FALSE)
    )
  }
  
  if(dist == "exp") {
    return(
      function(x) pexp(q = x,
                       rate = param.fit$estimate["rate"],
                       lower.tail = FALSE)
    )
  }
  
  if(dist == "lnorm") {
    return(
      function(x) plnorm(q = x,
                         meanlog = param.fit$estimate["meanlog"],
                         sdlog = param.fit$estimate["sdlog"],
                         lower.tail = FALSE)
    )
  }
  
  if(dist == "logis") {
    return(
      function(x) plogis(q = x,
                         location = param.fit$estimate["location"],
                         scale = param.fit$estimate["scale"],
                         lower.tail = FALSE)
    )
  }
}