# Tau-Process-Inference

### res.tau.hat_func()
This function calculates the estimated local Kendall's tau process, between a binary variable and a continuous variable subject to right censoring, up to a specified time point. Given a confidence level, this function produces the pointwise confidence interval of the local Kendall’s tau process evaluated at the specified time. The methods were proposed by Tai, Wang and Wells. <br>
(Tai, Y.C., Wang, W., & Wells, M.. (2022). Kendall's Tau for Two-Sample Inference Problems. https://arxiv.org/abs/2207.14445) <br>

#### Arguments
`X`: a non-empty numeric vector of group indicators, coded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of observed times <br>
`delta`: the censoring indicator coded as 0 if censored; as 1 if failed <br>
`t.star`: the specified time point (no greater than the minimum of the largest failure points in the two samples) <br>
`alpha`: type I error, default 0.05 <br>

#### Value
A list containing the following components <br>
`tau.hat`: the estimated value of the local Kendall’s &tau; process up to `t.star` <br>
`var.est`:  the variance of the local Kendall’s &tau;<br>
`ci.random`: the (1-alpha)*100% confidence interval of the local Kendall’s &tau; under random assignment <br>
`ci.fixed`: the (1-alpha)*100% confidence interval of the local Kendall’s &tau; under fixed assignment <br>
`p.est`: p-value using the estimated variance <br>
`p.0`: p-value using the null variance estimate <br>
`prop.compared`: the proportion of the concordance/discordance relationships being calculated up to `t.star` <br>

#### Example
The dataset was obtained from Pagano and Gauvreau (2000: exercise 9, page 512). There were 86 patients who received surgery to remove tumors. After the surgery, 48 patients were assigned to the placebo treatment (Group 0) and 38 patients were assigned to chemotherapy (Group 1). The variable under comparison is the time to first recurrence of tumor. The estimated local Kendall’s is plotted up to `t.star` = 59 (months). <br>

```
res.tau.hat_func(X = cancer$group,
                 observed.time = cancer$time,
                 delta = cancer$censor,
                 t.star = 30)
                 
$tau.hat
[1] 0.1316818

$var.est
[1] 0.01682179

$ci
[1] -0.1225232  0.3858868

$prop.compared
[1] 0.8321354
```


### res.tau.process_func()
This function wraps the procedure to calculate the local Kendall's &tau; process at the specified time point. 

#### Arguments
`X`: a non-empty numeric vector of group indicators, coded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of observed times <br>
`delta`: the censoring indicator coded as 0 if censored; as 1 if failed <br>
`t.star`: the vector of specified time points (no greater than the minimum of the largest failure points in the two samples) <br>

#### Value
It returns a vector of the estimated local Kendall's &tau; process at the specified time points. Furthermore, it automatically draws two plots, where the two Kaplan-Meier curves are on the left and the estimated tau process is on the right. 

#### Example
```
res.tau.process_func(X = cancer$group,
                     observed.time = cancer$time,
                     delta = cancer$censor,
                     t.star = quantile(cancer$time, seq(0, 1, by = 0.1)))

[1]  0.000000000 -0.059178744  0.042874396  0.001811594  0.024120606  0.105456868  0.125783386  0.114345016  0.131681805  0.145064589
[11]  0.145064589
```
![github_example_bladder](https://user-images.githubusercontent.com/9900943/180813073-e3430a0a-17a6-40ca-91a0-8c6354d7fb84.png)

### imputed.tau.hat_func2()
This function impute the concordane and discordance relationship within the provided interval. 

#### Arguments
`X`: a non-empty numeric vector of group indicators, coded as 0 or 1 <br>
`observed.time`: a non-empty numeric vector of observed times <br>
`delta`: the censoring indicator coded as 0 if censored; as 1 if failed <br>
`t.min`: The left endpoint of the interval being imputed <br>
`t.max`: The right endpoint of the interval being imputed <br>
`dist.0`: the assumed distribution of group 0. ("exp", "weibull", "lnorm", "logis") <br>
`dist.1`: the assumed distribution of group 1. ("exp", "weibull", "lnorm", "logis") <br>

#### Value
It returns the value of the imputed local Kendall's &tau; at `t.max`.

## Remark
The dependency packages include `survival` and `parmsurvfit`.

## Reference
1. Pagano, M., & Gauvreau, K. (2000). Principles of biostatistics. Australia: Duxbury. <br>
2. Tai, Y.C., Wang, W., & Wells, M.. (2022). Kendall's Tau for Two-Sample Inference Problems. https://arxiv.org/abs/2207.14445
