# Tau-Process-Inference

### res.tau.hat_func()
This function calculates the estimated local Kendall's tau process, between a binary variable and a continuous variable subject to right censoring, up to a specified time point. Given a confidence level, this function produces the pointwise confidence interval of the local Kendall’s tau process evaluated at the specified time. The methods were proposed by Tai, Wang and Wells. <br>

#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of observed failure times <br>
`delta`: the censoring indicator coded as 0 if censored; as 1 if failed <br>
`t.star`: the specified time point <br>
`alpha`: type I error, default 0.05 <br>

#### Value
A list containing the following components <br>
`tau.hat`: the estimated value of the local Kendall’s &tau; process up to `t.star` <br>
`var.est`:  the variance of the local Kendall’s &tau;<br>
`ci`: the (1-alpha)*100% confidence interval of the local Kendall’s &tau;<br>

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

### res.tau.process_func()
This function wrap the procedure to calculate the local Kendall's &tau; on th given vector of specified time points. 

#### Arguments
`X`: a non-empty numeric vector of group indicators, encoded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of observed failure times <br>
`delta`: the censoring indicator coded as 0 if censored; as 1 if failed <br>
`t.star`: the specified time point <br>

#### Value
The left hand side of the plot is the Kaplan-Meier curves of the controal and treatment groups. The right hand side is the line plot of the proposed tau process. 

#### Example

## Remark
The dependency packages include `survival`.
