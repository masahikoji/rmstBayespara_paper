---
title: 'rmstBayespara: An R Package for Bayesian Parametric Restricted Mean Survival
  Time'
tags:
- R
- restricted mean survival time
- Bayesian analysis
- Weibull distribution
- "log-normal distribution"
date: "28 Octber 2024"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Keisuke Hanada
  equal-contrib: true
  affiliation: 1
- name: Masahiro Kojima
  orcid: "0000-0003-0867-7692"
  equal-contrib: true
  affiliation: 2
bibliography: paper.bib
aas-doi: "10.3847/xxxxx <- update this with the DOI from AAS once you know it."
aas-journal: "The Journal of Open Source Software <- The name of the AAS journal."
affiliations:
- name: Osaka University, Japan
  index: 1
- name: The Institute of Statistical Mathematics, Japan
  index: 2
---

# Summary

The restricted mean survival time (RMST) has recently been widely used due to its straightforward interpretation as a mean measure. Until now, there has not been a package that derives the distribution of RMST using parametric Bayesian methods. The `rmstBayepara` package enables parametric Bayesian analysis considering covariates and heterogeneity. It provides posterior samples of RMST, shrinkage-adjusted RMST for each cluster, and RMST under specific covariate conditions. The summarized results include the mean, median, mode, and 95% credible interval for the mean. The analysis considering heterogeneity includes methods that add random effects to the linear combination of covariates and apply frailty to the hazard. The rmstBayepara package can analyze exponential, Weibull, log-normal, and log-logistic models. It also implements model selection using the Widely Applicable Information Criterion (WAIC).


# Statement of need

The restricted mean survival time (RMST) has recently been widely used due to its straightforward interpretation as a mean measure [@uno2014moving;@kim2017restricted]. The RMST is obtained by integrating the survival function up to a specified time limit. In non-parametric analysis methods, such as the Kaplan-Meier method, the survival function is estimated, and the RMST can be obtained by calculating the sum of the areas of each segment. In frequentist statistics, methods have been proposed for analyzing the influence of covariates and the heterogeneity between clusters[@andersen2003generalised;@tian2014predicting]. In Bayesian methods, non-parametric approaches for calculating RMST have been suggested by Zhang and Yin (2023) [@zhang2023bayesian,], but because they are non-parametric, it is not possible to perform analyses that adjust for covariates or consider heterogeneity.

For the RMST relating packages, `survRM2` provides RMST, the difference in RMST, the confidence interval for the difference, the ratio, and the confidence interval for the ratio. It also allows for covariate adjustment. The rmst function in the `RISCA` package can also calculate RMST. There is no a package of non-parametric Bayesian analysis.

The `rmstBayepara` package enables parametric Bayesian analysis considering covariates and heterogeneity. It provides posterior samples of RMST, shrinkage-adjusted RMST for each cluster, and RMST under specific covariate conditions. The summarized results include the mean, median, mode, and 95% credible interval for the mean. The analysis considering heterogeneity includes methods that add random effects to the linear combination of covariates and apply frailty to the hazard. The `rmstBayepara` package can analyze exponential, Weibull, log-normal, and log-logistic models. It also implements model selection using the Widely Applicable Information Criterion (WAIC) [@watanabe2010asymptotic].


<!--
# Results

We illustrate a real example of the Leukemia Survival Data in spBayesSurv, which include 1,043 patients of survival of acute myeloid leukemia. We focus on the "time", "cens", "age" and "district" variables.

```
data(LeukSurv)
dat<-LeukSurv
dat$status<-dat$cens
dat$arm[dat$age>=65]<-"Greater than or equal to 65 years old"
dat$arm[dat$age<65]<-"Less than 65 years old"
dat[c(1,1:5*200), c("time", "status", "arm", "age", "district")]
```

-->


<!--
```{r, ex-int, eval = knitr::is_html_output()}
data(LeukSurv)
dat<-LeukSurv
dat$status<-dat$cens
dat$arm[dat$age>=65]<-"Greater than or equal to 65 years old"
dat$arm[dat$age<65]<-"Less than 65 years old"
knitr::kable(dat[c(1,1:5*200), c("time", "status", "arm", "age", "district")], format="html",
             caption="Example dataset")
```
 -->
 

<!--
```{r, ex-static, eval = knitr::is_latex_output()}
data(LeukSurv)
dat<-LeukSurv
dat$status<-dat$cens
dat$arm[dat$age>=65]<-"Greater than or equal to 65 years old"
dat$arm[dat$age<65]<-"Less than 65 years old"
knitr::kable(head(dat[c(1,1:5*200), c("time", "status", "arm", "age", "district")]), format="latex",
             caption="Example dataset")
```
 -->



# Illustration of rmstBayespara

We introduce how to use `rmstBayespara`. The package `rmstBayespara` can be installed and loaded using the following code.
```
library(rmstBayespara)
library(spBayesSurv)
```

The following main functions are demonstrated in this article:

`brm_surv`: Inference on the Bayesian regression models for parametric survival time.

`rmstpara`: Restricted mean survival time via parametric models.

The function provides the inference of Bayesian regression models using Stan for parametric survival time.

**Usage**: The `brm_surv` function is called using following syntax.
```
brm_surv(time, cnsr, var, rvar, family = "exponential", random = "fixed", data,
  iter = 2000, warmup = 1000, seed = NA, chains = 4)
```

**Arguments**: The main arguments of `brm_surv` function is tabulated as below:

- `time`: name of time variable in data. Need to set character.
- `cnsr`: name of censor variable in data. Need to set character.
- `var`: vector of covariate names in data. Need to set character.
- `rvar`: name of random effect in data. Need to set character.
- `family`: A description of the response distribution and link function to be used in the model. 'exponential', 'Weibull', 'log-normal', and 'log-logistic' can be selected.
- `random`: A description of random effect. 'fixed', 'normal', and 'frailty' are available.
- `data`: An object of class data.frame (or one that can be coerced to that class) containing data of all variables used in the model.


**Value**: The `brm_surv` function returns a list of list of an object of class `brmsfit` or `stanfit` (see \CRANpkg{rstan} and \CRANpkg{brms}), sampling values from posterior distribution, leave-one-out cross-validation, and widely applicable information criterions.

- `fit`: An object of class `brmsfit` or `stanfit`.
- `post_sample`: Extract posterior samples of specified parameters.
- `waic`: Widely applicable information criterion.
- `loo`: Efficient approximate leave-one-out cross-validation.

The function provides the value of RMST via parametric models, such as exponential, Weibull, log-normal, and log-logistic models.

**Usage**: The `rmstpara` function is called using following syntax.
```
rmstpara(tau, var, rvar = NA, shape = NA, sigma = NA, 
         family = "exponential", random = "fixed")
```

**Arguments**: The main arguments of `rmstpara` function is tabulated as below:

- `tau`: A value of pre-specified evaluation time point.
- `var`: A vector of covariate values.
- `rvar`: A vector of frailty effects. It is necessary when log-normal frailty and log-logistic frailty models.
- `shape`: A vector of shape parameters. It is necessary when Weibull and log-logistic models.
- `sigma`: A vector of standard error parameters. It is necessary when log-normal model.
- `family`: A description of the response distribution and link function to be used in the model. 'exponential', 'Weibull', 'log-normal', and 'log-logistic' can be selected.
- `random`: A description of random effect. 'fixed', 'normal', and 'frailty' are available.

**Value**: The `rmstpara` function returns a value or vector of RMST via the specific parametric model.

## Example data
We illustrate a real example of the Leukemia Survival Data in `spBayesSurv`, which include 1,043 patients of survival of acute myeloid leukemia. We focus on the "time", "cens", "age" and "district" variables.
"District" can be regarded as a cluster, and the impact of mixed effects and frailty can be considered.

```
data(LeukSurv)
dat<-LeukSurv
dat$status<-dat$cens
dat$arm[dat$age>=65]<-"Greater than or equal to 65 years old"
dat$arm[dat$age<65]<-"Less than 65 years old"
dat[c(1,1:5*200), c("time", "status", "arm", "age", "district")]
```


## Sample code
First of all, we calculate the naive Weibull model with two variables `arm` and `sex`.
```
fit_x <- brm_surv(time="time", cnsr="1-status", 
                  var=c("factor(arm)", "factor(sex)"), 
                  data=dat, family="Weibull", random="fixed")
```


The WAIC and loo can be illustrated as below:

```               
fit_x$waic
fit_x$loo
```

||Estimate|SE|
| ---- | ---- | ---- |
|elpd_waic|-6056.5|73.6|
|p_waic|4.8|0.4|
|waic|12113.0|147.1|

||Estimate|SE|
| ---- | ---- | ---- |
|elpd_loo|-6056.5|73.6|
|p_loo|4.8|0.4|
|looic|12113.0|147.1|

Although main purpose of `rmstBayespara` is to present group differences in RMST, the Bayesian parametric model also can illustrate the posterior distribution of covariate parameters.

```
estimate <- apply(fit_x$post_sample, 2, summary)
round(estimate, 4)
boxplot(fit_x$post_sample, cex.axis=0.8)
```

||b_intercept|b_factor(arm)|b_factor(sex)|shape|
| ---- | ---- | ---- | ---- | ---- |
|Min.|5.4637|1.1233|-0.5284|0.5114|
|1st Qu.|5.7489|1.5012|-0.1710|0.5446|
|Median|5.8212|1.5867|-0.0865|0.5542|
|Mean|5.8219|1.5860|-0.0873|0.5547|
|3rd Qu. |5.8930|1.6692|-0.0058|0.5648|
|Max.|6.2475|2.0499|0.3431|0.6007|

The difference in RMST between arms of posterior samples from Weibull model can be obtained as below:

```
ps_x<-fit_x$post_sample
rmst_diff_x<-as.numeric(
  rmstpara(tau=100, var=ps_x[,'b_intercept']+ps_x[,'b_factor(arm)'], 
           shape=ps_x[,'shape'], family="Weibull",random="fixed") - 
    rmstpara(tau=100, var=ps_x[,'b_intercept'],
             shape=ps_x[,'shape'], family="Weibull",random="fixed")
  )
c(summary(rmst_diff_x),quantile(rmst_diff_x,c(0.025,0.975)))
hist(rmst_diff_x, breaks=100, freq=F)
```

|Min.|1st Qu.|Median|Mean|3rd Qu.|Max.|2.5%|97.5%|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|10.18390|13.89870|14.77579|14.77010|15.59717|19.53334|12.25126|17.27399|

The 95% confidence interval can be obtained from the 2.5% and 97.5% results.

Next, the analysis of mixed effects model can be implemented as follows.　By setting random="normal", it is possible to adapt to a mixed effects model.
```         
fit_x_r <- brm_surv(time="time", cnsr="1-status", var=c("factor(arm)", "factor(sex)"), 
                    rvar="district", data=dat, family="Weibull", random="normal")
fit_x_r$waic
fit_x_r$loo
```
||Estimate|SE|
| ---- | ---- | ---- |
|elpd_waic|-6050.0|73.8|
|p_waic|21.0|1.6|
|waic|12100.0|147.6|

||Estimate|SE|
| ---- | ---- | ---- |
|elpd_loo|-6050.1|73.8|
|p_loo|21.0|1.6|
|looic|12100.1|147.6|

The results of the difference in RMST with mixed effect, taking heterogeneity into account, can be calculated as follows.

```
ps_x_r<-fit_x_r$post_sample
rmst_diff_x_r<-as.numeric(
  rmstpara(tau=100, var=ps_x_r[,'b_intercept']+ps_x_r[,'b_factor(arm)'],
           rvar=ps_x_r[,'district'], shape=ps_x_r[,'shape'], 
           family="Weibull", random="normal") - 
    rmstpara(tau=100, var=ps_x_r[,'b_intercept'], 
           rvar=ps_x_r[,'district'], shape=ps_x_r[,'shape'], 
           family="Weibull", random="normal")
  )
c(summary(rmst_diff_x_r),quantile(rmst_diff_x_r,c(0.025,0.975)))
```

|Min.|1st Qu.|Median|Mean|3rd Qu.|Max.|2.5%|97.5%|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|8.837232|12.074202|12.917372|12.974916|13.844946|17.774116|10.443351|15.619410|

The results of the difference in RMST for distinct 1, taking heterogeneity into account, can be calculated as follows.
```
rmst_diff_x_r_1<-as.numeric(
  rmstpara(tau=100, var=ps_x_r[,'b_intercept']+ps_x_r[,'b_factor(arm)'], 
           shape=ps_x_r[,'shape'], rvar=ps_x_r[,'sd_district[1]'], 
           family="Weibull",random="normal") - 
    rmstpara(tau=100, var=ps_x_r[,'b_intercept'], 
             shape=ps_x_r[,'shape'], rvar=ps_x_r[,'sd_district[1]'], 
             family="Weibull",random="normal")
  )
c(summary(rmst_diff_x_r_1),quantile(rmst_diff_x_r_1,c(0.025,0.975)))
```

|Min.|1st Qu.|Median|Mean|3rd Qu.|Max.|2.5%|97.5%|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|0.2712211|0.4375881|0.4887243|0.4945920|0.5460238|0.9143291|0.3492974|0.6731644|

Finally, frailty effects model can be implemented as follows.
```         
fit_x_f <- brm_surv(time="time", cnsr="1-status", var=c("factor(arm)", "factor(sex)"), 
                    rvar="district", data=dat, family="Weibull", random="frailty")
fit_x_f$waic
fit_x_f$loo
```
||Estimate|SE|
| ---- | ---- | ---- |
|elpd_waic|-6050.0|73.9|
|p_waic|22.0|1.7|
|waic|12100.0|147.7|

||Estimate|SE|
| ---- | ---- | ---- |
|elpd_loo|-6050.1|73.9|
|p_loo|22.1|1.7|
|looic|12100.2|147.7|

The results of the difference in RMST with frailty effect, taking heterogeneity into account, can be calculated as follows.

```
ps_x_f<-fit_x_f$post_sample
rmst_diff_x_f<-as.numeric(
  rmstpara(tau=100, var=ps_x_f[,'b_intercept']+ps_x_f[,'b_factor(arm)'],
           rvar=ps_x_f[,'sd_district'], shape=ps_x_f[,'shape'], 
           family="Weibull", random="frailty") - 
    rmstpara(tau=100, var=ps_x_f[,'b_intercept'], 
           rvar=ps_x_f[,'sd_district'], shape=ps_x_f[,'shape'], 
           family="Weibull", random="frailty")
  )
c(summary(rmst_diff_x_f),quantile(rmst_diff_x_f,c(0.025,0.975)))
```

|Min.|1st Qu.|Median|Mean|3rd Qu.|Max.|2.5%|97.5%|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|0.08846989|0.68157877|1.02491134|1.13388377|1.44250681|4.85410016|0.27837340|2.70644121|

The results of the difference in RMST for distinct 1, taking heterogeneity into account, can be calculated as follows.
```
rmst_diff_x_f_1<-as.numeric(
  rmstpara(tau=100, var=ps_x_f[,'b_intercept']+ps_x_f[,'b_factor(arm)'], 
           shape=ps_x_f[,'shape'],rvar=ps_x_f[,"sd_district[1]"],  
           family="Weibull",random="frailty") - 
    rmstpara(tau=100, var=ps_x_f[,'b_intercept'], 
             shape=ps_x_f[,'shape'], rvar=ps_x_f[,"sd_district[1]"], 
             family="Weibull",random="frailty")
  )
c(summary(rmst_diff_x_f_1),quantile(rmst_diff_x_f_1,c(0.025,0.975)))
```
|Min.|1st Qu.|Median|Mean|3rd Qu.|Max.|2.5%|97.5%|
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
|-25.039281591|8.385407514|10.359507614|9.357048164|11.463067134|14.833317067|0.005389485|12.930702706|

# Acknowledgements

# References
