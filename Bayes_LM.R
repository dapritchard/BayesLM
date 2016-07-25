
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
#' @param decomp_method One of either "chol" or "eigen", specifying whether the
#'     multivariate normal sampling function should use the Cholesky
#'     decomposition or the eigen decomposition.


bayes_lm_r <- function(n                = 100L,
                       p                = 15L,
                       nsamp            = 1e4,
                       prop_nonzero     = 0.2,
                       true_beta_sd     = 2,
                       true_sigma       = 2,
                       sd_x             = 2,
                       print_stats      = FALSE,
                       time_sections    = TRUE,
                       write_samples    = FALSE,
                       samples_file_loc = "Samples_R.dat",
                       decomp_method    = "chol") {

    # Ensure that user-supplied arguments are of the right type and have legal
    # values.  source check_valid_input() if necessary.
    if (! exists("check_valid_input")) {
        source("Check_Valid_Input.R")
    }
    check_valid_input(n, p, nsamp, prop_nonzero, true_beta_sd, true_sigma, sd_x,
                      print_stats, time_sections, write_samples,
                      samples_file_loc, decomp_method)


    # Set values of model objects ----------------------------

    # Sample values of beta.  Want to evenly space (100 * prop_nonzero)% of the
    # values from beta to be nonzero.  Indices calculated by evenly spacing each
    # nonzero value on the real line, and then truncating to an integer.
    # n_nonzero + 1 is the number of spaces surrounding the nonzero values.
    n_nonzero <- as.integer(prop_nonzero * p)
    nonzero_idx <- as.integer((p / (n_nonzero + 1)) * 1:n_nonzero) + 1L
    true_beta <- replace(rep(0, p),
                         nonzero_idx,
                         rnorm(length(nonzero_idx), sd=true_beta_sd))

    # Sample data
    X <- matrix(rnorm(n * p, 0, sd_x), nrow=n, ncol=p)
    y <- X %*% true_beta + rnorm(n, 0, true_sigma)

    # Set the priors; see Hoff pgs. 154-155 for the meanings of the priors.
    # Note: we are implicitely specifying the mean hyperparameter for beta to be
    # 0 by ommitting the term in the Gibbs sampler conditional mean calculation.
    Sigma_inv_0 <- diag(p)
    nu_0 <- 1
    sigma_sq_0 <- 1

    # Write true values of beta, sigma^{-2} to the first row of output file
    if (write_samples) {
        write(c(true_beta, true_sigma^(-2)), samples_file_loc, p + 1L)
    }


    # Preliminary calculations  ------------------------------

    tXX <- crossprod(X)
    tXy <- crossprod(X, y)
    nu_sigma_sq_0 <- nu_0 * sigma_sq_0
    shapeval <- (nu_0 + n) / 2


    # Sampler loop -------------------------------------------

    # Conditionally allocate memory to store samples in
    if (print_stats) {
        out_beta <- matrix(0, p, nsamp)
        out_gamma <- vector("numeric", nsamp)
    }

    # Initial value for gamma
    gamma <- 1

    # Set pointer to specified multivariate normal sampling method
    samp_mvnorm <- switch(decomp_method,
                          chol  = quote( drop( mvtnorm::rmvnorm(1L, m, V, "chol") ) ),
                          eigen = quote( MASS::mvrnorm(1L, m, V) ),
                          stop("illegal decomp_method argument", call.=FALSE))

    # Clock timer start
    elapsed_inv <- 0
    elapsed_samp <- 0
    start_time <- proc.time()

    for (s in 1:nsamp) {

        # Sample beta.  case: don't time inverse / sampling sections
        if (!time_sections) {
            V <- solve(Sigma_inv_0 + (tXX * gamma))
            m <- V %*% tXy * gamma
            beta <- eval(samp_mvnorm)
        }
        # case: do time inverse / sampling sections
        else {
            elapsed_inv <- elapsed_inv +
                system.time(V <- solve(Sigma_inv_0 + (tXX * gamma)), FALSE)[3L]
            m <- V %*% tXy * gamma
            elapsed_samp <- elapsed_samp +
                system.time(beta <- eval(samp_mvnorm), FALSE)[3L]
        }

        # Sample gamma
        SSR <- crossprod(y - (X %*% beta))
        rateval <- (nu_sigma_sq_0 + SSR) / 2
        gamma <- rgamma(1L, shapeval, rateval)

        # Conditionally store data in memory / write to file
        if (print_stats) {
            out_beta[, s] <- beta
            out_gamma[s] <- gamma
        }
        if (write_samples) {
            write(c(beta, gamma), samples_file_loc, p + 1L, TRUE)
        }
    }

    # Calculate elapsed time
    elapsed_overall <- (proc.time() - start_time)[3L]


    # Print summary statistics -------------------------------

    if (print_stats) {

        # beta table of true values and posterior quantiles
        beta_quants <- apply(out_beta, 1, quantile, c(0.025, 0.5, 0.975))
        beta_tab <- data.frame(true_beta,
                               t(beta_quants))
        names(beta_tab) <- c("true", "2.5%", "50%", "97.5%")

        # gamma table of true values and posterior quantiles
        gamma_quants <- quantile(out_gamma, c(0.025, 0.5, 0.975))
        gamma_tab <- data.frame(1 / true_sigma^2,
                                matrix(gamma_quants, nrow=1))
        names(gamma_tab) <- c("true", "2.5%", "50%", "97.5%")

        # Write beta and gamma tables to console
        cat("\n",
            "true beta and posterior quantiles:\n",
            "----------------------------------\n",
            "\n", sep="")
        print(beta_tab, digits=4, row.names=FALSE)

        cat("\n",
            "true gamma and posterior quantiles:\n",
            "-----------------------------------\n",
            "\n", sep="")
        print(gamma_tab, digits=4, row.names=FALSE)
        cat("\n")
    }

    c(inverse  = elapsed_inv,
      sampling = elapsed_samp,
      overall  = elapsed_overall)
}
