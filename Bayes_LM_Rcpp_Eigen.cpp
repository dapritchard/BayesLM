
// [[Rcpp::depends(RcppEigen)]]

#include <fstream>
#include <RcppEigen.h>
#include <ctime>               // clock_gettime

#include "Stats_Fcns_Eigen.h"  // sample_beta, mvrnorm_chol, mvrnorm_eigen


#define NANO_MULT 0.000000001  // nano multiplier, i.e. 1e-9

#define OVERALL   0  // index for overall time elapsed
#define INVERSE   1  // index for time elapsed calculating matrix inverse
#define SAMP_NORM 2  // index for time elapsed sampling from normal distribution


#define CLOCK_START(idx) clock_gettime(CLOCK_MONOTONIC, &start[idx]);
#define CLOCK_STOP(idx) do {						\
	clock_gettime(CLOCK_MONOTONIC, &finish[idx]);			\
	elapsed[idx] += (finish[idx].tv_sec - start[idx].tv_sec);	\
	elapsed[idx] += (finish[idx].tv_nsec - start[idx].tv_nsec) * NANO_MULT; \
    } while (0)


// Comparator function for use in qsort
int compare_dbl(const void* a, const void* b) {
    return (* (double*) a) - (* (double*) b) < 0 ? -1 : 1;
}



// [[Rcpp::export]]
Rcpp::NumericVector bayes_lm_rcpp_eigen_(int n,
					int p,
					int nsamp,
					double prop_nonzero,
					double true_beta_sd,
					double true_sigma,
					double sd_x,
					bool print_stats,
					bool write_samples,
					Rcpp::CharacterVector samples_file_loc,
					char decomp_method) {


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
	samples_file.open(samples_file_loc[0].begin());
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
    double* bcurr;         // beta current data location
    double* bend;          // 1 past beta data end location

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
	samples_file.open(samples_file_loc[0].begin(), std::fstream::app);
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
	SSR = err.squaredNorm();
	scaleval = 2 / (nu_sigma_sq_0 + SSR);
	gamma = R::rgamma(shapeval, scaleval);

	// Conditionally store data in memory / write to file
	if (write_samples) {
	    bcurr = beta.data();
	    bend = bcurr + p;
	    for (; bcurr < bend; bcurr++) {
		samples_file << *bcurr << " ";
	    }
	    samples_file << gamma << "\n";
	}
	if (print_stats) {
	    out_beta.col(s) = beta;
	    out_gamma(s) = gamma;
	}

	// Allow user to interrupt and return to the R REPL
	if (s % 1000 == 0) {
	    Rcpp::checkUserInterrupt();
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
	
	Rcpp::Rcout << "\n"
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
		    << "write_samples:  " << write_samples << "\n"
		    << "samples_file_loc:  " << samples_file_loc << "\n"
		    << "decomp_method:  " << decomp_method << "\n";

	// Set printing of fields to be a fixed format with precision 4
	Rcpp::Rcout.setf(std::ios::fixed, std::ios::floatfield);
	Rcpp::Rcout.precision(4);

	Rcpp::Rcout << "\n"
		    << "Elapsed time:\n"
		    << "-------------\n"
		    << "Inverse:          " << elapsed[INVERSE] << "\n"
		    << "Sampling normal:  " << elapsed[SAMP_NORM] << "\n"
		    << "Overall:          " << elapsed[OVERALL] << "\n"
		    << "\n";

	// Set printing of matrices
	Eigen::IOFormat matprint(4, 0, "  ", "\n", "  ", "", "", "");
	Eigen::IOFormat gamprint(4, 0, "   ", "\n", "   ", "", "", "");
	
	Rcpp::Rcout << "true beta     2.5%     50%     97.5%\n"
		    << "------------------------------------\n"
		    << table_beta_quant.format(matprint) << "\n"
		    << "\n"
		    << " true gam     2.5%     50%     97.5%\n"
		    << "------------------------------------\n"
		    << table_gamma_quant.format(gamprint) << "\n"
		    << "\n";
    }
  
    // Return computational time 
    return Rcpp::NumericVector::create(Rcpp::Named("inverse")  = elapsed[INVERSE],
    				       Rcpp::Named("mvnsamp") = elapsed[SAMP_NORM],
    				       Rcpp::Named("overall")  = elapsed[OVERALL]);
}
