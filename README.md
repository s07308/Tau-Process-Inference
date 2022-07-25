# Tau-Process-Inference

### res.tau.hat_func()


#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of data <br>
`delta`: the status indicator. Typically, 0: censored, 1: failed <br>
`t.star`: the specified cut point <br>
`alpha`: type I error, default 0.05 <br>

#### Value
A list containing the following components <br>
`tau.hat`: the estimated value of &tau;<sub>b</sub> <br>
`var.est`: the variance of the estimator of &tau;<sub>b</sub><br>
`ci`: the (1-alpha)*100% confidence interval of &tau;<sub>b</sub><br>

#### Example
The dataset is obtained from the book "Survival Analysis Techniques for Censored and Truncated Data" (Klein, John P., Moeschberger, Melvin L., 2003). It recoreded the time to infection for patients receiving Kidney Dialysis. <br>

```
time.knot <- quantile(KD$time, seq(0, 1, by = 0.1))
tau.process <- numeric(length = length(time.knot))

for(k in seq_along(time.knot)) {
  tau.process.fit <- res.tau.hat_func(X = KD$treatment,
                                      observed.time = KD$time,
                                      delta = KD$delta,
                                      t.star = time.knot[k])
  tau.process[k] <- tau.process.fit$tau.hat
}

plot(time.knot, tau.process, type = "b")
```
![github_example_kidney](https://user-images.githubusercontent.com/9900943/180728285-496d57d8-044e-4fd1-aca3-5c808185bd39.png)

## Remark
The dependency packages include `survival`.
