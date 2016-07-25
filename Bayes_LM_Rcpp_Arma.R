
#' Bayesian linear model Gibbs sampler; see Hoff pages 154-155 for mathematical
#' details.
#'
#' The main purpose of the program is as an exercise to approximate the
#' posterior distribution of the Bayesian linear model using a Gibbs sampler
#' with an R implementation, and then compare the speed of the implementation
#' with that of a comparable implementation written in C++.
#'
#' In more detail, we are estimating p(\beta, \gamma | X, y) using a Gibbs
#' sampler for the model y = X\beta + \epsilon.  The values X, \beta, and y are
#' sampled, and then the posterior distributions of \beta and \gamma are
#' approximated through the Gibbs sampler.  Summary statistics may be written to
#' stdout, and the MCMC samples may be written to file.
#'
#' The model hyperparameters are specified with fixed values as follows:
#'
#'     beta_0:     zero vector
#'     sigma_0:    identity matrix
#'     nu_0:       1
#'     sigma_0^2:  1
#'
#'
#' @param n A numeric value no less than 1 specifying the number of observations
#'     in the model.  Non-integer values are truncated to integer values.
#'
#' @param p A numeric value no less than 1 specifying the number of predictor
#'     variables in the model.  Non-integer values are truncated to integer
#'     values.
#'
#' @param nsamp A numeric value no less than 1 specifying the number of scans to
#'     perform in the Gibbs sampler.  Non-integer values are truncated to
#'     integer values.
#'
#' @param prop_nonzero A numeric value in the range (0, 1] specifying the
#'     proportion of nonzero predictor coefficients in the model.
#'
#' @param true_beta_sd A positive numeric value; each nonzero predictor
#'     coefficient in the true model is independently sampled from from a N(0,
#'     true_beta_sd^2) distribution.
#'
#' @param true_sigma A positive numeric value; the model outcome vector y is
#'     sampled from a conditional distribution y | X\beta ~ N(0, true_sigma^2),
#'     where X is the matrix of predictor variables, and \beta is a vector of
#'     variable coefficients.
#'
#' @param sd_x A positive numeric value; the model predictor variables are
#'     independently sampled from a N(0, sd_x^2) distribution.
#'
#' @param print_stats One of either TRUE or FALSE, specifying whether a printout
#'     of the true (sampled) predictor coefficient values and approximations of
#'     the 2.5%, 50%, and 97.5% quantile levels is written to the console.
#'
#' @param write_samples One of either TRUE or FALSE, specifying whether the
#'     samples generated from the Gibbs sampler should be written to file.  The
#'     first row provides the values for the true beta vector and inverse of the
#'     true sigma^2, and the following rows provide the samples.  The location
#'     of the file is specified by samples_file_loc.
#'
#' @param samples_file_loc A character string specifying the location of the
#'     file to which the samples generated from the Gibbs sampler should be
#'     written (ignored if write_samples is FALSE).
#'
#' @param method One of either "chol" or "eigen", specifying whether the
#'     multivariate normal sampling function should use the Cholesky
#'     decomposition or the eigen decomposition.


bayes_lm_rcpp_arma <- function(n                = 100L,
                               p                = 15L,
                               nsamp            = 1e4,
                               prop_nonzero     = 0.2,
                               true_beta_sd     = 2,
                               true_sigma       = 2,
                               sd_x             = 2,
                               print_stats      = FALSE,
                               write_samples    = FALSE,
                               samples_file_loc = "Samples_Rcpp_Arma.dat",
                               decomp_method    = "chol") {

    # Ensure that user-supplied arguments are of the right type and have legal
    # values.  source check_valid_input() if necessary.
    if (! exists("check_valid_input")) {
        source("Check_Valid_Input.R")
    }
    # We don't have an option for time_sections in this function - as opposed to
    # 'bayes_lm_r' - so just pass TRUE to 'check_valid_input'
    check_valid_input(n, p, nsamp, prop_nonzero, true_beta_sd, true_sigma, sd_x,
                      print_stats, TRUE, write_samples, samples_file_loc,
                      decomp_method)

    # To enable casting to a char in callee function
    decomp_method <- switch(decomp_method, "chol" = "c", "eigen" = "e")

    # Compile source code and load into memory
    if (! exists("bayes_lm_rcpp_arma_")) {
        Rcpp::sourceCpp("Bayes_LM_Rcpp_Arma.cpp")
    }

    # Call the underlying workhorse function
    bayes_lm_rcpp_arma_(n, p, nsamp, prop_nonzero, true_beta_sd, true_sigma, sd_x,
                        print_stats, write_samples, samples_file_loc, decomp_method)
}
