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
The dataset is obtained from Pagano and Gauvreau (2000: exercise 9, page 512). There were 86 patients who received surgery to remove tumors. After the surgery, 48 patients were assigned to the placebo treatment (Group 0) and 38 patients were assigned to chemotherapy (Group 1). The variable under comparison is the time to first recurrence of tumor. The estimated local Kendall’s is plotted up to `t.star` = 59 (months). <br>

```

```

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
```

```
## Remark
The dependency packages include `survival`.
