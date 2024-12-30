# rmstBayespara

Bayesian regression models using 'Stan' for restricted mean survival time. The package implement the model estimation described in Hanada and Kojima (2024)

Hanada, K., & Kojima, M. (2024). Bayesian Parametric Methods for Deriving Distribution of Restricted Mean Survival Time. arXiv e-prints, arXiv-2406 (Submitted).

## brm_surv

A function "brm_surv" is Bayesian regression models using stan for parametric survival time. Exponential, Weibull, log-normal, and log-logistic model with fixed-effect, random-effect and frailty-effect can be available.

## rmstpara

A function "rmstpara" is calculating restricted mean survival time via parametric models. Exponential, Weibull, log-normal and log-logistic models are available.


## Example
Given times, statuses, arm, a baseline characteristic (sex), and a cluster (district), Bayesian regression models and the restricted mean survival time can be calculated as follows:

```R
d <- data.frame(time=1:100,
      status=sample(0:1, size=100, replace=TRUE),
      arm=sample(c("t", "c"), size=100, replace=TRUE),
      sex=sample(1:2, size=100, replace=TRUE),
      district=sample(1:5, size=100, replace=TRUE)
     )
head(d)
fit_x_r <- brm_surv(time="time", cnsr="1-status",
                    var=c("factor(arm)", "factor(sex)"),
                    rvar="district", data=d,
                    family="Weibull", random="frailty")
fit_x_r$fit
fit_x_r$post_sample
fit_x_r$waic
fit_x_r$loo
ps_x_r<-fit_x_r$post_sample
rmst_x_r<-rmstpara(tau=100, var=ps_x_r[,"b_intercept"]+ps_x_r[,"b_factor(arm)"],
                   shape=ps_x_r[,"shape"], rvar=ps_x_r[,"sd_district"],
                   family="Weibull",random="frailty")
rmst_x_r

```
