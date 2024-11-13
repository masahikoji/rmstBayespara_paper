#'Bayesian regression models using 'Stan' for survival time
#' @description A function of Bayesian regression models using stan for parametric survival time. Exponential, Weibull, log-normal, and log-logistic model with fixed-effect, random-effect and frailty-effect can be available.
#'
#' @name brm_surv
#'
#' @param time name of time variable in data. Need to set character.
#' @param cnsr name of censor variable in data. Need to set character.
#' @param var vector of covariate names in data. Need to set character.
#' @param rvar name of random effect in data. Need to set character.
#' @param family A description of the response distribution and link function to be used in the model. 'exponential', 'Weibull', 'log-normal', and 'log-logistic' can be selected.
#' @param random A description of random effect. 'fixed', 'normal', and 'frailty' are available.
#' @param data An object of class data.frame (or one that can be coerced to that class) containing data of all variables used in the model.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup (aka burnin) iterations. This also specifies the number of iterations used for stepsize adaptation, so warmup draws should not be used for inference. The number of warmup should not be larger than iter and the default is iter/2.
#' @param seed The seed for random number generation to make results reproducible. If NA (the default), Stan will set the seed randomly.
#' @param chains Number of Markov chains (defaults to 4).
#'
#' @return A list of an object of class brmsfit or stanfit (see rstan and brms), sampling values from posterior distribution, leave-one-out cross-validation, and widely applicable information criterions.
#'
#' @examplesIf interactive()
#' d <- data.frame(time=1:100,
#'       status=sample(0:1, size=100, replace=TRUE),
#'       arm=sample(c("t", "c"), size=100, replace=TRUE),
#'       sex=sample(1:2, size=100, replace=TRUE),
#'       district=sample(1:5, size=100, replace=TRUE)
#'     )
#' head(d)
#' fit_x_r <- brm_surv(time="time", cnsr="1-status",
#'                     var=c("factor(arm)", "factor(sex)"),
#'                     rvar="district", data=d,
#'                     family="Weibull", random="frailty"
#'                     )
#' fit_x_r$fit
#' fit_x_r$post_sample
#' fit_x_r$waic
#' fit_x_r$loo
#'
#'
#' @export
brm_surv <- function(time, cnsr, var, rvar, family="exponential", random="fixed", data, iter=2000, warmup=1000, seed=NA, chains=4){
  v <- base::paste(var, collapse = "+")
  if(random=="fixed"){
    fval <- base::paste(time, "|", "cens(", cnsr, ") ~ ", v, sep="")
    f <- stats::as.formula(fval)
  }else if(random=="normal" | random=="frailty"){
    fval <- base::paste(time, "|", "cens(", cnsr, ") ~ ", v, "+ (1|", rvar, ")", sep="")
    f <- stats::as.formula(fval)
  }else{
    stop("'random' variable must be set to 'fixed', 'normal', or 'frailty'.")
  }
  sdat <- brms::make_standata(f, data)

  if(family=="exponential"){
    if(random=="fixed" | random=="normal"){
      x <- brms::brm(f,
                     data=data,
                     family = exponential(),
                     iter=iter,
                     warmup=warmup,
                     seed=seed,
                     chains=chains)
    }else if(random=="frailty"){
      exp_frail <- "
        // generated with brms 2.21.0
        functions {
        }
        data {
          int<lower=1> N;  // total number of observations
          vector[N] Y;  // response variable
          array[N] int<lower=-1,upper=2> cens;  // indicates censoring
          int<lower=1> K;  // number of population-level effects
          matrix[N, K] X;  // population-level design matrix
          int<lower=1> Kc;  // number of population-level effects after centering
          // data for group-level effects of ID 1
          int<lower=1> N_1;  // number of grouping levels
          int<lower=1> M_1;  // number of coefficients per level
          array[N] int<lower=1> J_1;  // grouping indicator per observation
          // group-level predictor values
          vector[N] Z_1_1;
          int prior_only;  // should the likelihood be ignored?
        }
        transformed data {
          matrix[N, Kc] Xc;  // centered version of X without an intercept
          vector[Kc] means_X;  // column means of X before centering
          for (i in 2:K) {
            means_X[i - 1] = mean(X[, i]);
            Xc[, i - 1] = X[, i] - means_X[i - 1];
          }
        }
        parameters {
          vector[Kc] b;  // regression coefficients
          real Intercept;  // temporary intercept for centered predictors
          vector<lower=0>[M_1] sd_1;  // group-level standard deviations
          vector<lower=0>[N_1] v;
        }
        transformed parameters {
          real lprior = 0;  // prior contributions to the log posterior
          lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        }
        model {
          v~gamma(1/sd_1[1],1/sd_1[1]);
          // likelihood including constants
          if (!prior_only) {
            // initialize linear predictor term
            vector[N] mu = rep_vector(0.0, N);
            mu += Intercept + Xc * b;
            mu = exp(mu);
            for (n in 1:N) {
              // add more terms to the linear predictor
              mu[n] = mu[n]*v[J_1[n]] ;
            }
            for (n in 1:N) {
            // special treatment of censored data
              if (cens[n] == 0) {
                target += exponential_lpdf(Y[n] | inv(mu[n]));
              } else if (cens[n] == 1) {
                target += exponential_lccdf(Y[n] | inv(mu[n]));
              } else if (cens[n] == -1) {
                target += exponential_lcdf(Y[n] | inv(mu[n]));
              }
            }
          }
          // priors including constants
          target += lprior;
        }
        generated quantities {
          // actual population-level intercept
          real b_Intercept = Intercept - dot_product(means_X, b);
          vector[N] log_lik;
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          mu = exp(mu);
          for (n in 1:N) {
            // add more terms to the linear predictor
            mu[n] = mu[n]*v[J_1[n]] ;
          }
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
              log_lik[n] = exponential_lpdf(Y[n] | inv(mu[n]));
            } else if (cens[n] == 1) {
              log_lik[n] = exponential_lccdf(Y[n] | inv(mu[n]));
            } else if (cens[n] == -1) {
              log_lik[n] = exponential_lcdf(Y[n] | inv(mu[n]));
            }
          }
        }
        "
      base::message(crayon::red("Compiling Stan program...\n"))
      exp_frail_model <- rstan::stan_model(model_code = exp_frail)
      base::message(crayon::red("Start sampling \n"))
      x <- rstan::sampling(exp_frail_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }
  }else if(family=="Weibull"){
    if(random=="fixed" | random=="normal"){
      x <- brms::brm(f,
                     data=data,
                     family = weibull(),
                     iter=iter,
                     warmup=warmup,
                     seed=seed,
                     chains=chains)
    }else if(random=="frailty"){
      weibull_frail <- "
      // generated with brms 2.21.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        array[N] int<lower=-1,upper=2> cens;  // indicates censoring
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=0> shape;  // shape parameter
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        vector<lower=0>[N_1] v;
      }
      transformed parameters {
        real lprior = 0;  // prior contributions to the log posterior
        lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        lprior += gamma_lpdf(shape | 0.01, 0.01);
      }
      model {
        v~gamma(1/sd_1[1],1/sd_1[1]);
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          mu = exp(mu);
          for (n in 1:N) {
            // add more terms to the linear predictor
            mu[n] = mu[n]/pow(v[J_1[n]],1/shape);
          }
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
              target += weibull_lpdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            } else if (cens[n] == 1) {
              target += weibull_lccdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            } else if (cens[n] == -1) {
              target += weibull_lcdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            }
          }
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] log_lik;
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b;
        mu = exp(mu);
        for (n in 1:N) {
          // add more terms to the linear predictor
          mu[n] = mu[n]/pow(v[J_1[n]],1/shape);
        }
        for (n in 1:N) {
        // special treatment of censored data
            if (cens[n] == 0) {
              log_lik[n] = weibull_lpdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            } else if (cens[n] == 1) {
              log_lik[n] = weibull_lccdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            } else if (cens[n] == -1) {
              log_lik[n] = weibull_lcdf(Y[n] | shape, mu[n] / tgamma(1 + 1 / shape));
            }
        }
      }
      "
      weibull_frail_model <- rstan::stan_model(model_code = weibull_frail)
      x <- rstan::sampling(weibull_frail_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }
  }else if(family=="log-normal"){
    if(random=="fixed" | random=="normal"){
      x <- brms::brm(f,
                     data=data,
                     family = lognormal(),
                     iter=iter,
                     warmup=warmup,
                     seed=seed,
                     chains=chains)
    }else if(random=="frailty"){
      ln_frail <- "
      // generated with brms 2.21.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        array[N] int<lower=-1,upper=2> cens;  // indicates censoring
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=0> sigma;  // dispersion parameter
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        vector<lower=0>[N_1] v;
      }
      transformed parameters {
        real lprior = 0;  // prior contributions to the log posterior
        lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        lprior += student_t_lpdf(sigma | 3, 0, 2.5)
          - 1 * student_t_lccdf(0 | 3, 0, 2.5);
      }
      model {
        v~gamma(1/sd_1[1],1/sd_1[1]);
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
      //        target += lognormal_lpdf(Y[n] | mu[n], sigma)+log(v[J_1[n]])+(v[J_1[n]]-1)*lognormal_lccdf(Y[n] | mu[n], sigma);
              target += lognormal_lpdf(Y[n] | mu[n], sigma)+log(v[J_1[n]])+(v[J_1[n]]-1)*log(1-lognormal_cdf(Y[n] | mu[n], sigma));
            } else if (cens[n] == 1) {
              target += v[J_1[n]]*log(1-lognormal_cdf(Y[n] | mu[n], sigma));
      //        target += v[J_1[n]]*lognormal_lccdf(Y[n] | mu[n], sigma);
            }
          }
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] log_lik;
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b;
        for (n in 1:N) {
        // special treatment of censored data
            if (cens[n] == 0) {
              log_lik[n] = lognormal_lpdf(Y[n] | mu[n], sigma)+log(v[J_1[n]])+(v[J_1[n]]-1)*lognormal_lccdf(Y[n] | mu[n], sigma);
            } else if (cens[n] == 1) {
              log_lik[n] = v[J_1[n]]*lognormal_lccdf(Y[n] | mu[n], sigma);
            }
        }
      }
      "
      base::message(crayon::red("Compiling Stan program...\n"))
      ln_frail_model <- rstan::stan_model(model_code = ln_frail)
      base::message(crayon::red("Start sampling \n"))
      x <- rstan::sampling(ln_frail_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }
  }else if(family=="log-logistic"){
    if(random=="fixed"){
      ll_fixed <- "
      // generated with brms 2.21.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        array[N] int<lower=-1,upper=2> cens;  // indicates censoring
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=1> shape;  // shape parameter
      }
      transformed parameters {
        real lprior = 0;  // prior contributions to the log posterior
        lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        lprior += gamma_lpdf(shape | 0.01, 0.01);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          mu=exp(mu);
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
              target += loglogistic_lpdf(Y[n] | mu[n], shape);
            } else if (cens[n] == 1) {
              target += log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            }
          }
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] log_lik;
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b;
        mu=exp(mu);
        for (n in 1:N) {
        // special treatment of censored data
          if (cens[n] == 0) {
            log_lik[n] = loglogistic_lpdf(Y[n] | mu[n], shape);
          } else if (cens[n] == 1) {
            log_lik[n] = log(1-loglogistic_cdf(Y[n] | mu[n], shape));
          }
        }
      }
      "
      base::message(crayon::red("Compiling Stan program...\n"))
      ll_fixed_model <- rstan::stan_model(model_code = ll_fixed)
      base::message(crayon::red("Start sampling \n"))
      x <- rstan::sampling(ll_fixed_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }else if(random=="normal"){
      ll_normal <- "
      // generated with brms 2.21.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        array[N] int<lower=-1,upper=2> cens;  // indicates censoring
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=1> shape;  // shape parameter
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        array[M_1] vector[N_1] z_1;  // standardized group-level effects
      }
      transformed parameters {
        vector[N_1] v;  // actual group-level effects
        real lprior = 0;  // prior contributions to the log posterior
        v = (sd_1[1] * (z_1[1]));
        lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        lprior += gamma_lpdf(shape | 0.01, 0.01);
        lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
          - 1 * student_t_lccdf(0 | 3, 0, 2.5);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          for (n in 1:N) {
            // add more terms to the linear predictor
            mu[n] += v[J_1[n]] * Z_1_1[n];
          }
          mu = exp(mu);
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
              target += loglogistic_lpdf(Y[n] | mu[n], shape);
            } else if (cens[n] == 1) {
              target += log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            }
          }
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] log_lik;
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b;
        for (n in 1:N) {
          // add more terms to the linear predictor
          mu[n] += v[J_1[n]] * Z_1_1[n];
        }
        mu=exp(mu);
        for (n in 1:N) {
        // special treatment of censored data
          if (cens[n] == 0) {
            log_lik[n] = loglogistic_lpdf(Y[n] | mu[n], shape);
          } else if (cens[n] == 1) {
            log_lik[n] = log(1-loglogistic_cdf(Y[n] | mu[n], shape));
          }
        }
      }
      "
      base::message(crayon::red("Compiling Stan program...\n"))
      ll_normal_model <- rstan::stan_model(model_code = ll_normal)
      base::message(crayon::red("Start sampling \n"))
      x <- rstan::sampling(ll_normal_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }else if(random=="frailty"){
      ll_frail <- "
      // generated with brms 2.21.0
      functions {
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        array[N] int<lower=-1,upper=2> cens;  // indicates censoring
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=1> shape;  // shape parameter
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        vector<lower=0>[N_1] v;
      }
      transformed parameters {
        real lprior = 0;  // prior contributions to the log posterior
        lprior += student_t_lpdf(Intercept | 3, 5.3, 2.5);
        lprior += gamma_lpdf(shape | 0.01, 0.01);
      }
      model {
        v~gamma(1/sd_1[1],1/sd_1[1]);
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          mu = exp(mu);
          for (n in 1:N) {
            // add more terms to the linear predictor
            mu[n] = mu[n]*v[J_1[n]] ;
          }
          for (n in 1:N) {
          // special treatment of censored data
            if (cens[n] == 0) {
              target += loglogistic_lpdf(Y[n] | mu[n], shape)+log(v[J_1[n]])+(v[J_1[n]]-1)*log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            } else if (cens[n] == 1) {
              target += v[J_1[n]]*log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            }
          }
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
        vector[N] log_lik;
        vector[N] mu = rep_vector(0.0, N);
        mu += Intercept + Xc * b;
        mu = exp(mu);
        for (n in 1:N) {
        // special treatment of censored data
            if (cens[n] == 0) {
              log_lik[n] = loglogistic_lpdf(Y[n] | mu[n], shape)+log(v[J_1[n]])+(v[J_1[n]]-1)*log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            } else if (cens[n] == 1) {
              log_lik[n] = v[J_1[n]]*log(1-loglogistic_cdf(Y[n] | mu[n], shape));
            }
        }
      }
      "
      base::message(crayon::red("Compiling Stan program...\n"))
      ll_frail_model <- rstan::stan_model(model_code = ll_frail)
      base::message(crayon::red("Start sampling \n"))
      x <- rstan::sampling(ll_frail_model,
                           data=sdat,
                           iter=iter,
                           warmup=warmup,
                           seed=seed,
                           chains=chains)
    }
  }else{
    stop("'family' variable must be set to 'exponential', 'Weibull', 'log-normal', 'log-logistic'.")
  }

  post_sample <- brms::as_draws_matrix(x)
  cname <- colnames(post_sample)
  n_var <- 1 + length(var) # intercept + covariate
  cname[1:n_var] <- paste("b_", c("intercept", var), sep="")
  if(family!="exponential"){
    n_var <- n_var + 1 # shape or sigma
  }
  if(random=="normal" | random=="frailty"){
    n_var <- n_var + 1
    if(brms::is.brmsfit(x)){
      cname[n_var-1] <- base::paste("sd_", rvar, sep="")
    }else{
      cname[n_var] <- base::paste("sd_", rvar, sep="")
    }
    n_r <- base::length(base::unique(data[,rvar]))
    cname[n_var+1:n_r] <- base::paste("sd_", rvar, "[", base::sort(base::unique(data[,rvar])), "]", sep="")
  }else{
    n_r <- 0
  }
  post_sample <- post_sample[,1:(n_var+n_r)]
  if(!brms::is.brmsfit(x)){
    if(random=="fixed"){
      post_sample <- post_sample[,1:n_var]
    }else{
      post_sample <- post_sample[,c(length(var)+1, 1:length(var), (length(var)+2):n_var, n_var+1:n_r)]
    }
  }
  colnames(post_sample) <- cname[1:(n_var+n_r)]

  x_loo <- loo::loo(x)
  if(brms::is.brmsfit(x)){
    x_waic <- brms::waic(x)
  }else{
    x_waic <- loo::waic(loo::extract_log_lik(x))
  }

  return(list(fit=x, post_sample=post_sample, loo=x_loo, waic=x_waic))
}
