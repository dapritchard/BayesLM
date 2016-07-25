
// case: when linking against Rmath as a stand-alone library
#ifdef MATHLIB_STANDALONE

#include <Eigen/Dense>
#include <Rmath.h>      // rnorm
#include <stdexcept>    // runtime_error

/* Creates the same alias when using Rmath as a stand-alone library as when we
 * are calling this code through R via Rcpp.  This allows us to use the same
 * function call for both scenarios.
 */
namespace R {
inline double rnorm(double mu, double sigma) { return ::rnorm(mu, sigma); }
}

/* case: if not linking against Rmath as a stand-alone library, then assumes
 * that we are linking against R via Rcpp
 */
#else
#include <RcppEigen.h>
#endif


#define TOLERANCE 0.000000000001  // 1e-12




// Statistical functions -------------------------------------------------------

/* Generates values for beta where the number of nonzero values in beta is given
 * by floor(prop_nonzero * p) and are evenly spaced along the vector.  The
 * nonzero entries of beta are sampled from a N(0, true_beta_sd) distribution.
 */

Eigen::VectorXd sample_beta(int p, double prop_nonzero, double true_beta_sd) {

    Eigen::VectorXd beta(p);  // Allocate memory to store beta vector
    int n_nonzero;            // Number of nonzero elements of beta
    double stepsize;          // Fractional amount of zero entries between nonzero
    double next_idx_dbl;      // Fractional index of next nonzero entry
    int next_idx_int;         // Integer index of next nonzero entry

    /* Note that in the following casts from double to int we need not consider
     * overflow since p is known to be a representable number, and both
     * prop_nonzero * p and next_idx_dbl are no bigger than p.
     */
    n_nonzero = (int) (prop_nonzero * p);
    stepsize = p / (double) (n_nonzero + 1);
    next_idx_dbl = stepsize;
    next_idx_int = (int) next_idx_dbl;

    /* Each iteration sets the value of beta(j).  Add a machine error smudge
     * factor mainly to prevent the last element from being erroneously
     * considered a nonzero element due to the sum of stepsize being less than p
     * - 1 due to machine error.
     */
    for (int j = 0; j < p; j++) {

	/* case: a nonzero value.  Sample beta(j).  Next, calculate where the
	 * next value would be placed if we could evenly specify it on the real
	 * line, and then set the next index for a nonzero value to be the floor
	 * of the true evenly spaced value.
	 */
	if (j == next_idx_int) {
	    beta(j) = R::rnorm(0, true_beta_sd);
	    next_idx_dbl += stepsize;
	    next_idx_int = (int) (next_idx_dbl + TOLERANCE);
	}
	// case: a zero value for beta
	else {
	    beta(j) = 0;
	}
    }

    return beta;
}




/* Returns an nrow by ncol matrix such that the elements in the matrix are
 * randomly sampled from independent N(0, sd^2) distributions.
 */

Eigen::MatrixXd matr_randn(int nrow, int ncol, double sd) {

    Eigen::MatrixXd X(nrow, ncol);
    double* curr = X.data();
    double* end = curr + X.size();

    for ( ; curr < end; curr++) {
	*curr = R::rnorm(0, sd);
    }

    return X;
}




/* Generates a random sample from a N(m, V) distribution using an eigen
 * decomposition of V.  Note that the function uses the fact the the eigenvalues
 * are provided by eig_sym in increasing order.
 *
 * PRE: assumes that V is symmetric and positive-definite, and that m has the
 * same number of elements as V has rows.
 */

Eigen::VectorXd mvrnorm_eigen(Eigen::VectorXd& m, Eigen::MatrixXd& V) {

    int p = m.size();
    Eigen::VectorXd eigval;       // Store eigenvalues of V
    Eigen::VectorXd root_eigval;  // Store square root of e-values of V
    Eigen::MatrixXd eigvec;       // Store eigenvectors of V
    Eigen::VectorXd mvnsamp;      // Store a random sample from indep N(0, 1)'s

    // Perform eigenvalue decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_decomp(V);
    if (eig_decomp.info() != Eigen::Success) {
    	throw std::runtime_error("Error in eigen decomposition");
    }
    eigval = eig_decomp.eigenvalues();
    eigvec = eig_decomp.eigenvectors();

    /* Calculate elementwise square root of eigenvalues.  Ensure beforehand that
     * eigenvalues that are negative due to machine error are set to zero.
     */
    double* curr = eigval.data();
    double* end = curr + p;
    // Iterate through eigenvalues until we reach positive eigenvalues
    for ( ; (curr != end) && (*curr < 0); curr++) {
	*curr = 0;
    }
    root_eigval = eigval.array().sqrt();

    /* Eigen decomposition provides the relation V = PDP^{T} where P is an
     * orthogonal matrix where the columns are eigenvectors and D is a diagonal
     * matrix composed of the corresponding eigenvalues.  It follows that for
     * independent normal r.v. z:
     *
     *     Var[P D^{1/2} z] = P D^{1/2} Var[z] D^{1/2} P^{T}
     *                      = P D^{1/2} I D^{1/2} P^{T}
     *                      = PDP^{-1}
     *                      = V
     */
    return m + (eigvec * root_eigval.cwiseProduct(matr_randn(p, 1, 1)));
}




/* Generates a random sample from a N(m, V) distribution using an Cholesky
 * decomposition of V.
 *
 * PRE: assumes that V is symmetric and positive-definite, and that m has the
 * same number of elements as V has rows.
 */

Eigen::VectorXd mvrnorm_chol(Eigen::VectorXd& m, Eigen::MatrixXd& V) {

    // Perform Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> R(V);
    if (R.info() != Eigen::Success) {
    	throw std::runtime_error("Error in Cholesky decomposition");
    }

    /* Cholesky decomposition provides the relation V = RR^{T} so that for
     * independent normal r.v. z:
     *
     *     Var[Rz] = R Var[z] R^{T} = RIR^{T} = RR^{T} = V
     */
    return m + (R.matrixL() * matr_randn(m.size(), 1, 1));
}




/* Calculate the prob-th empirical quantile for each column of X and return as
 * elements of a Eigen::VectorXd.  This is an implementation of the R function
 * stats::quantile under the default choice of empirical quantile.  See the help
 * page for stats::quantile for more details (the default there is type 7).  The
 * formula is (for n the number of observations in the data):
 *
 *     F^{-1}(prob) = ((1 - gamma) * x_j) + (gamma * x_{j+1})
 *
 * where
 *
 *     j     = floor(prob * (n - 1))
 *     gamma = (prob * (n - 1)) - floor(prob * (n - 1))
 *
 * PRE: assumes X is sorted within each column.
 */

Eigen::VectorXd quantile(Eigen::MatrixXd& X, double prob) {

    if ((prob <= 0) || (prob >= 1)) {
	throw std::runtime_error("prob must be in the range (0, 1)");
    }
    else if (X.rows() < 2) {
	// Need 2+ observations to prevent formula from invalid index call
	throw std::runtime_error("X must have 2 or more rows of data");
    }

    double desired_idx;  // where you might choose point if you had a continuous
			 // set of points
    int below_idx;       // the closest point below desired point
    int above_idx;       // the closest point above desired point
    double gamma;        // linear interpolation value

    Eigen::VectorXd quant_vals(X.cols());  // container for quantile values

    desired_idx = prob * (X.rows() - 1);
    below_idx   = (int) desired_idx;
    above_idx   = below_idx + 1;
    gamma       = desired_idx - below_idx;

    // Each iteration calculates the prob-th empirical quantile for the j-th col
    unsigned int ncols = X.cols();
    for (unsigned int j = 0; j < ncols; j++) {
	quant_vals(j) = ((1 - gamma) * X(below_idx, j)) + (gamma * X(above_idx, j));
    }

    return quant_vals;
}




/* Overloaded version of Eigen::VectorXd quantile(Eigen::MatrixXd& X, double
 * prob).  See comments for this function for details of how the empirical
 * quantile is defined.
 */

double quantile(Eigen::VectorXd& X, double prob) {

    if ((prob <= 0) || (prob >= 1)) {
	throw std::runtime_error("prob must be in the range (0, 1)");
    }
    else if (X.rows() < 2) {
	// Need 2+ observations to prevent formula from invalid index call
	throw std::runtime_error("X must have length no less than 2");
    }

    double desired_idx;  // where you might choose point if you had a continuous
			 // set of points
    int below_idx;       // the closest point below desired point
    int above_idx;       // the closest point above desired point
    double gamma;        // linear interpolation value

    desired_idx = prob * (X.rows() - 1);
    below_idx   = (int) desired_idx;
    above_idx   = below_idx + 1;
    gamma       = desired_idx - below_idx;

    return ((1 - gamma) * X(below_idx)) + (gamma * X(above_idx));
}
