
/* Bayesian linear model Gibbs sampler using Eigen library; see Hoff pages
 * 154-155 for mathematical details.
 *
 * The main purpose of the program is as an exercise to use the Amadillo library
 * in an MCMC sampler setting, and also to test the speed of the sampler under
 * various options.
 *
 * In more detail, we are estimating p(\beta, \gamma | X, y) using a Gibbs
 * sampler for the model y = X\beta + \epsilon.  The values X, \beta, and y are
 * sampled, and then the posterior distributions of \beta and \gamma are
 * approximated through the Gibbs sampler.  Summary statistics may be written to
 * stdout, and the MCMC samples may be written to file.
 *
 * The model hyperparameters are specified with fixed values as follows:
 *
 *     beta_0:     zero vector
 *     sigma_0:    identity matrix
 *     nu_0:       1
 *     sigma_0^2:  1
 */

/* The following command line arguments can be provided.  The desired value of
 * the argument should follow the argument specifier, either immediately after
 * or separated by whitespace.  
 *
 * @param -n An integer value no less than 1 specifying the number of
 *     observations in the model.  Non-integer values are truncated to integer
 *     values.
 *
 * @param -p A numeric value no less than 1 specifying the number of predictor
 *     variables in the model.  Non-integer values are truncated to integer
 *     values.
 *
 * @param -nsamp A numeric value no less than 1 specifying the number of scans
 *     to perform in the Gibbs sampler.  Non-integer values are truncated to
 *     integer values.
 *
 * @param -prop_nonzero A numeric value in the range (0, 1] specifying the
 *      proportion of nonzero predictor coefficients in the model.
 *
 * @param -true_beta_sd A positive numeric value; each nonzero predictor
 *     coefficient in the true model is independently sampled from from a N(0,
 *     true_beta_sd^2) distribution.
 *
 * @param -true_sigma A positive numeric value; the model outcome vector y is
 *     sampled from a conditional distribution y | X\beta ~ N(0, true_sigma^2),
 *     where X is the matrix of predictor variables, and \beta is a vector of
 *     variable coefficients.
 *
 * @param -sd_x A positive numeric value; the model predictor variables are
 *     independently sampled from a N(0, sd_x^2) distribution.
 *
 * @param -print_stats One of either true or false, specifying whether a
 *     printout of the true (sampled) predictor coefficient values and
 *     approximations of the 2.5%, 50%, and 97.5% quantile levels is written to
 *     the console.
 *
 * @param -write_samples One of either true or false, specifying whether the
 *     samples generated from the Gibbs sampler should be written to file.  The
 *     first row provides the values for the true beta vector and inverse of the
 *     true sigma^2, and the following rows provide the samples.  The location
 *     of the file is specified by samples_file_loc.
 *
 * @param -samples_file_loc A character string specifying the location of the
 *     file to which the samples generated from the Gibbs sampler should be
 *     written (ignored if write_samples is false).
 *
 * @param -write_ctime One of either true or false, specifying whether the
 *     computation time taken to perform calculate a matrix inverse, sample from
 *     the normal distribution, and the overall time should by written to
 *     file. The location of the file is specified by ctime_file_loc.
 *
 * @param -ctime_file_loc A character string specifying the location of the file
 *     to which the computational time should be written (ignored if write_ctime
 *     is false).
 *
 * @param -decomp_method One of either "chol" or "eigen", specifying whether the
 *     multivariate normal sampling function should use the Cholesky
 *     decomposition or the eigen decomposition.
 *
 * @param -seed A nonnegative integer specifying a value that the RNG is to be
 *     seeded with.
 */

/* Compile using e.g.
 *
 *        g++ Bayes_LM_Eigen.cpp Parse_Args.cpp Stats_Fcns_Eigen.cpp  \
 *            -DMATHLIB_STANDALONE -Wall -g3 -O3 -lRmath              \
 *            -o bayes_lm_eigen
 *
 * This requires having the Eigen headers and Rmath library available.
 */

/* CAVEAT: this file may require a POSIX system to compile due to the use of
 * clock_gettime() for the clock timer.  Users of other systems may require a
 * modification of that part of the code.
 */


#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Rmath.h>      // rgamma
#include <ctime>        // time
#include <stdexcept>    // runtime_error
#include <cstdlib>      // qsort

#include "Parse_Args.h"        // parse_args, FILENAME_MAXLEN
#include "Stats_Fcns_Eigen.h"  // sample_beta, mvrnorm_chol, mvrnorm_eigen


#define TIME_MULTIPLIER  79    // arbitrary prime number to obtain second seed num
#define NANO_MULT 0.000000001  // nano multiplier, i.e. 1e-9

#define OVERALL   0  // index for overall time elapsed
#define INVERSE   1  // index for time elapsed calculating matrix inverse
#define SAMP_NORM 2  // index for time elapsed sampling from normal distribution


// Track time of inverse, sampling from normal, and overall
#ifndef NO_TIMER
#define CLOCK_START(idx) clock_gettime(CLOCK_MONOTONIC, &start[idx]);
#define CLOCK_STOP(idx) do {						\
	clock_gettime(CLOCK_MONOTONIC, &finish[idx]);			\
	elapsed[idx] += (finish[idx].tv_sec - start[idx].tv_sec);	\
	elapsed[idx] += (finish[idx].tv_nsec - start[idx].tv_nsec) * NANO_MULT; \
    } while (0)
#else
#define CLOCK_START(idx) (void) start[idx];  // prevent unused variable warning
#define CLOCK_STOP(idx) (void) finish[idx];  // prevent unused variable warning
#endif


// Comparator function for use in qsort
int compare_dbl(const void* a, const void* b) {
    return (* (double*) a) - (* (double*) b) < 0 ? -1 : 1;
}




// Default parameter specifications --------------------------------------------
 
// Specify size of data
int n = 100;  // number of observations
int p =  15;  // number of predictor variables

// Specify proportion of nonzero elements in beta
double prop_nonzero = 0.2;

// Specify beta coefficient standard deviation
double true_beta_sd = 2;

// Specify sd of eps_i in y_i = t(X_i) * beta + eps_i
double true_sigma = 2;

// Specify st. dev. in sampling predictors coeffs from indep N(0, sd_x^2) dists
double sd_x = 2;

// Specify number of MCMC scans
int nsamp = 1e4;

// Print sample median statistics
bool print_stats = false;

// Write computational time to output file
bool write_ctime = false;

// Computational time output file location
char ctime_file_loc[FILENAME_MAXLEN] = "Comp_Time_Eigen.dat";

// Write samples to file
bool write_samples = false;

// MCMC samples output file location
char samples_file_loc[FILENAME_MAXLEN] = "Samples_Eigen.dat";

// Specify whether to use Cholesky or eigen decomposition for sampling normals
char decomp_method = 'c';

// Specify seed for RNG
unsigned int seed = time(NULL);




//  Begin main -----------------------------------------------------------------

int main(int argc, char* argv[]) {

    // Read command-line arguments
    parse_args(argc, argv);

    // Set seed for random draws
    set_seed(seed, seed * TIME_MULTIPLIER);


    // Declare model data structures -----------------------

    Eigen::VectorXd true_beta;    // true coefficient vector

    Eigen::VectorXd y;            // response values
    Eigen::MatrixXd X;            // predictor coefficient matrix

    Eigen::VectorXd beta;         // current value of beta sample
    double gamma;                 // current value of gamma sample

    Eigen::MatrixXd Sigma_inv_0;  // inverse of beta variance hyperparam
    double nu_0;                  // hyperparam 1 for inverse-gamma prior
    double sigma_sq_0;            // hyperparam 2 for inverse-gamma prior

    Eigen::MatrixXd out_beta;     // memory for beta samples
    Eigen::VectorXd out_gamma;    // memory for gamma samples


    // Declare storage and timer data structures -----------

    std::ofstream ctime_file;         // sampler loop computational time file
    std::ofstream samples_file;       // samples file

    struct timespec start[3];         // store event starting time information
    struct timespec finish[3];        // store event ending time information
    double elapsed[3] = { 0, 0, 0 };  // tracks event cumulative elapsed time

    double* curr;                     // pointer steps through current val
    double* end;                      // pointer to mark one past the last val


    // Set values of model objects -------------------------

    // Set values of beta
    true_beta = sample_beta(p, prop_nonzero, true_beta_sd);

    // Sample data
    X = matr_randn(n, p, sd_x);
    y = (X * true_beta) + matr_randn(n, 1, true_sigma);

    /* Set the priors; see Hoff pgs. 154-155 for the meanings of the priors.
     * Note: we are implicitely specifying the mean hyperparameter for beta to
     * be 0 by ommitting the term in the Gibbs sampler conditional mean
     * calculation.
     */
    Sigma_inv_0.setIdentity(p, p);
    nu_0 = 1;
    sigma_sq_0 = 1;


    // Write param vals to file ----------------------------

    // Write true values of beta, sigma^{-2} to the first row of output file
    if (write_samples) {
	samples_file.open(samples_file_loc);
	curr = true_beta.data();
	end = curr + p;
	for ( ; (curr != end); curr++) {
	    samples_file << *curr << " ";
	}
	samples_file << 1 / (true_sigma * true_sigma) << "\n";
	samples_file.close();
    }


    // Preliminary calculations ----------------------------

    Eigen::MatrixXd tXX;   // value of X^{T} X
    Eigen::VectorXd tXy;   // value of X^{T} y
    double shapeval;       // shape parameter for gamma distribution samples
    double nu_sigma_sq_0;  // product of nu_0 and sigma^2_0

    tXX           = X.transpose() * X;
    tXy           = X.transpose() * y;
    nu_sigma_sq_0 = nu_0 * sigma_sq_0;
    shapeval      = (nu_0 + n) / 2;


    // Sampler object initialization -----------------------

    Eigen::MatrixXd V;     // variance of current beta sample
    Eigen::VectorXd m;     // mean of current beta sample
    Eigen::VectorXd err;   // model error, i.e. y - X \beta
    Eigen::MatrixXd iden;  // storing matrix identity
    double SSR;            // SSR (sum of squared errors)
    double scaleval;       // scale parameter for gamma distribution samples

    // Set pointer to desired multivariate normal sampling function
    Eigen::VectorXd (*samp_mvnorm)(Eigen::VectorXd&, Eigen::MatrixXd&);
    switch (decomp_method) {
    case 'c':
    	samp_mvnorm = &mvrnorm_chol;
    	break;
    case 'e':
    	samp_mvnorm = &mvrnorm_eigen;
    	break;
    default:
    	throw std::runtime_error("Illegal value of decomp_method");
    }

    // Conditionally allocate memory for samples
    if (print_stats) {
	out_beta.resize(p, nsamp);
	out_gamma.resize(nsamp, 1);
    }
    
    // Conditionally open samples file stream
    if (write_samples) {
	samples_file.open(samples_file_loc, std::fstream::app);
    }

    // Initial value for gamma
    gamma = 1;

    // Initialize identity matrix
    iden.setIdentity(p, p);


    // Sampler loop ----------------------------------------

    // Clock timer objects and initialization; requires POSIX system
    CLOCK_START(OVERALL);

    for (int s = 0; s < nsamp; s++) {
	
	// Sample beta
	CLOCK_START(INVERSE)
	    // Leverage the fact that we have a p.d. matrix to obtain inverse
	    V = (Sigma_inv_0 + (gamma * tXX)).llt().solve(iden);
	CLOCK_STOP(INVERSE);
	m = gamma * V * tXy;
	CLOCK_START(SAMP_NORM)
	    beta = samp_mvnorm(m, V);
	CLOCK_STOP(SAMP_NORM);

	// Sample gamma
	err = y - (X * beta);
	// root_SSR = norm(err);
	SSR = err.dot(err);
	scaleval = 2 / (nu_sigma_sq_0 + SSR);
	gamma = rgamma(shapeval, scaleval);

	// Conditionally store data in memory / write to file
	if (write_samples) {
	    curr = beta.data();
	    end = curr + p;
	    for (; curr < end; curr++) {
		samples_file << *curr << " ";
	    }
	    samples_file << gamma << "\n";
	}
	if (print_stats) {
	    out_beta.col(s) = beta;
	    out_gamma(s) = gamma;
	}
    }

    // Calculate elapsed time
    CLOCK_STOP(OVERALL);


    // Print summary statistics ----------------------------

    if (print_stats) {

	// Allocate memory for tables with cols true values and quantiles
	Eigen::MatrixXd table_beta_quant(p, 4);
	Eigen::MatrixXd table_gamma_quant(1, 4);

	// Sort data in preperation for quantile function
	out_beta.transposeInPlace();
	double* col_start = out_beta.data();
	for (int j = 0; j < p; j++) {
	    qsort(col_start, nsamp, sizeof(double), compare_dbl);
	    col_start += nsamp;
	}
	qsort(out_gamma.data(), nsamp, sizeof(double), compare_dbl);

	// Fill in table entries for beta samples quantiles
	table_beta_quant.col(0) = true_beta;
	table_beta_quant.col(1) = quantile(out_beta, 0.025);
	table_beta_quant.col(2) = quantile(out_beta, 0.500);
	table_beta_quant.col(3) = quantile(out_beta, 0.975);

	// Fill in table entries for gamma samples quantiles
	table_gamma_quant(0, 0) = 1 / (true_sigma * true_sigma);
	table_gamma_quant(0, 1) = quantile(out_gamma, 0.025);
	table_gamma_quant(0, 2) = quantile(out_gamma, 0.500);
	table_gamma_quant(0, 3) = quantile(out_gamma, 0.975);
	
	std::cout << "\n"
		  << "Parameter specifications:\n"
		  << "-------------------------\n"
		  << "n:  " << n << "\n"
		  << "p:  " << p << "\n"
		  << "prop_nonzero:  " << prop_nonzero << "\n"
		  << "true_beta_sd:  " << true_beta_sd << "\n"
		  << "true_sigma:  " << true_sigma << "\n"
		  << "sd_x:  " << sd_x << "\n"
		  << "nsamp:  " << nsamp << "\n"
		  << "print_stats:  " << print_stats << "\n"
		  << "write_ctime:  " << write_ctime << "\n"
		  << "ctime_file_loc:  " << ctime_file_loc << "\n"
		  << "write_samples:  " << write_samples << "\n"
		  << "samples_file_loc:  " << samples_file_loc << "\n"
		  << "decomp_method:  " << decomp_method << "\n"
		  << "seed:  " << seed << "\n";

	// Set printing of fields to be a fixed format with precision 4
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(4);

	std::cout << "\n"
		  << "Elapsed time:\n"
		  << "-------------\n"
		  << "Inverse:          " << elapsed[INVERSE] << "\n"
		  << "Sampling normal:  " << elapsed[SAMP_NORM] << "\n"
		  << "Overall:          " << elapsed[OVERALL] << "\n"
		  << "\n";

	// Set printing of matrices
	Eigen::IOFormat matprint(4, 0, " ", "\n", " ", "", "", "");
	Eigen::IOFormat gamprint(4, 0, "   ", "\n", "   ", "", "", "");
	
	std::cout << "true beta     2.5%     50%     97.5%\n"
		  << "------------------------------------\n"
		  << table_beta_quant.format(matprint) << "\n"
		  << "\n"
		  << " true gam     2.5%     50%     97.5%\n"
		  << "------------------------------------\n"
		  << table_gamma_quant.format(gamprint) << "\n"
		  << "\n";
    }


    // Write computational time to output ------------------

    if (write_ctime) {
	ctime_file.open(ctime_file_loc, std::fstream::app);
	ctime_file << elapsed[INVERSE] << " "
		   << elapsed[SAMP_NORM] << " "
		   << elapsed[OVERALL] << "\n";
    }
  
    return 0;
}
