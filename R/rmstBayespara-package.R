#' The 'rmstBayespara' package.
#'
#' @description The parametric Bayes regression models using 'Stan' for restricted mean survival time (RMST).
#'            The package implement the model estimation described in Hanada and Kojima (2024).
#'            The posterior sample, Widely Applicable Information Criterion (WAIC), and
#'            Efficient Leave-One-Out Cross-Validation (loo) for some parametric models of
#'            survival data can be available. The RMST under these parametric models can be
#'            computed from the obtained posterior samples.
#'
#' @name rmstBayespara-package
#' @aliases rmstBayespara
#' @import methods
#' @import brms
#' @importFrom rstan stan_model
#' @importFrom rstan sampling
#' @importFrom loo loo
#' @importFrom loo waic
#' @importFrom loo extract_log_lik
#' @importFrom stats as.formula
#' @importFrom stats pnorm
#' @importFrom zipfR Ibeta
#' @importFrom zipfR Igamma
#'
#' @export brm_surv
#' @export rmstpara
#'
#' @references
#' Hanada, K., & Kojima, M. (2024). Bayesian Parametric Methods for Deriving Distribution of Restricted Mean Survival Time. arXiv preprint arXiv:2406.06071.
#'
"_PACKAGE"
