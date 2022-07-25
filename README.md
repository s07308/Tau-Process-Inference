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
