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


## Example data
We illustrate a real example of the Leukemia Survival Data in `spBayesSurv`, which include 1,043 patients of survival of acute myeloid leukemia. We focus on the "time", "cens", "age" and "district" variables.

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

Although main purpose of `rmstBayespara` is to present group differences in RMST, the Bayesian parametric model also can illustrate the posterior distribution of covariate parameters.

```
estimate <- apply(fit_x$post_sample, 2, summary)
round(estimate, 4)
boxplot(fit_x$post_sample, cex.axis=0.8)
```


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



Next, the analysis of mixed effects model can be implemented as follows.
```         
fit_x_r <- brm_surv(time="time", cnsr="1-status", var=c("factor(arm)", "factor(sex)"), 
                    rvar="district", data=dat, family="Weibull", random="normal")
fit_x_r$waic
fit_x_r$loo
```


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



Finally, frailty effects model can be implemented as follows.
```         
fit_x_f <- brm_surv(time="time", cnsr="1-status", var=c("factor(arm)", "factor(sex)"), 
                    rvar="district", data=dat, family="Weibull", random="frailty")
fit_x_f$waic
fit_x_f$loo
```



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




# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"


# Acknowledgements


# References