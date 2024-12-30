d <- data.frame(time=1:100,
     status=sample(0:1, size=100, replace=TRUE),
     arm=sample(c("t", "c"), size=100, replace=TRUE),
     sex=sample(1:2, size=100, replace=TRUE),
     district=sample(1:5, size=100, replace=TRUE)
   )

test_that(desc="Error1: family", {
  expect_error(brm_surv(time="time", cnsr="1-status",
                        var=c("factor(arm)", "factor(sex)"),
                        rvar="district", data=d,
                        family="weibull", random="frailty"),
               "'family' variable must be set to 'exponential', 'Weibull', 'log-normal', 'log-logistic'.")
})

test_that(desc="Error1: random", {
  expect_error(brm_surv(time="time", cnsr="1-status",
                        var=c("factor(arm)", "factor(sex)"),
                        rvar="district", data=d,
                        family="Weibull", random="abc"),
               "'random' variable must be set to 'fixed', 'normal', or 'frailty'.")
})

test_that(desc="Error2: family", {
  expect_error(rmstpara(tau=100, var=d$time, rvar=d$district, family="weibull", random="frailty"),
               "'family' variable must be set to 'exponential', 'Weibull', 'log-normal', 'log-logistic'.")
})

test_that(desc="Error2: shape", {
  expect_error(rmstpara(tau=100, var=d$time, rvar=d$district, family="Weibull", random="abc"),
               "'shape' variable need to calculate RMST when 'family' is 'Weibull'.")
})

test_that(desc="Error1: random", {
  expect_error(rmstpara(tau=100, var=d$time, shape=d$sex, rvar=d$district, family="Weibull", random="abc"),
               "'random' variable must be set to 'fixed', 'normal', or 'frailty'.")
})

