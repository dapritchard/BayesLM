
// case: when linking against Rmath as a stand-alone library
#ifdef MATHLIB_STANDALONE

#include <armadillo>
#include <R.h>        // rPsort
#include <Rmath.h>    // rnorm

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
#include <RcppArmadillo.h>
#endif


#define TOLERANCE 0.000000000001  // 1e-12




// Statistical functions -------------------------------------------------------

/* Generates values for beta where the number of nonzero values in beta is given
 * by floor(prop_nonzero * p) and are evenly spaced along the vector.  The
 * nonzero entries of beta are sampled from a N(0, true_beta_sd) distribution.
 */

arma::vec sample_beta(int p, double prop_nonzero, double true_beta_sd) {

    arma::vec beta(p);    // Allocate memory to store beta vector
    int n_nonzero;        // Number of nonzero elements of beta
    double stepsize;      // Fractional amount of zero entries between nonzero
    double next_idx_dbl;  // Fractional index of next nonzero entry
    int next_idx_int;     // Integer index of next nonzero entry

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
     * considered a nonzero element due to the sum of stepsize being less than
     * (p - 1) due to machine error.
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




/* Generates a random sample from a N(m, V) distribution using an eigen
 * decomposition of V.  Note that the function uses the fact the the eigenvalues
 * are provided by eig_sym in increasing order.
 *
 * PRE: assumes that V is symmetric and positive-definite, and that m has the
 * same number of elements as V has rows.
 */

arma::vec mvrnorm_eigen(arma::vec& m, arma::mat& V) {

    int p = m.n_elem;

    arma::vec eigval(p);       // Store eigenvalues of V
    arma::vec root_eigval(p);  // Store square root of e-values of V
    arma::mat eigvec(p, p);    // Store eigenvectors of V
    arma::vec mvnsamp(p);      // Store a random sample from indep N(0, 1)'s

    // Perform eigenvalue decomposition
    if (! eig_sym(eigval, eigvec, V)) {
	throw std::runtime_error("V is not positive definite\n");
    }

    /* Calculate elementwise square root of eigenvalues.  Ensure beforehand that
     * eigenvalues that are negative due to machine error are set to zero.
     */
    arma::vec::iterator eval_curr = eigval.begin(), eval_end = eigval.end();
    // Iterate through eigenvalues until we reach positive eigenvalues
    for ( ; (eval_curr != eval_end) && (*eval_curr < 0); eval_curr++) {
	*eval_curr = 0;
    }    
    root_eigval = sqrt(eigval);

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
    return m + (eigvec * (root_eigval % mvnsamp.randn()));
}




/* Generates a random sample from a N(m, V) distribution using an Cholesky
 * decomposition of V.
 *
 * PRE: assumes that V is symmetric and positive-definite, and that m has the
 * same number of elements as V has rows.
 */

arma::vec mvrnorm_chol(arma::vec& m, arma::mat& V) {

    int p = m.n_elem;
    arma::mat R(p, p);     // store Cholesky decomposition of V
    arma::vec mvnsamp(p);  // store a random sample from indep N(0, 1)'s

    /* Perform Cholesky decomposition.  Note: an exception is thrown if
     * decomposition fails.
     */
    R = chol(V, "lower");

    /* Cholesky decomposition with Armadillo option "lower" provides the
     * relation V = RR^{T} so that for independent normal r.v. z:
     *
     *     Var[Rz] = R Var[z] R^{T} = RIR^{T} = RR^{T} = V
     */
    return m + (R * mvnsamp.randn());
}




/* Calculates the empirical quantiles of the data pointed to by values and with
 * quantile levels corresponding to the indices provided by quant_below_idx(k)
 * and quant_below_idx(k) + 1 and interpolated by gamma(k) - see below for
 * definitions of quant_below_idx (corresponds to j) and gamma.  Return object
 * is the vector of quantiles.
 *
 * PRE: assumes that data_len >= 2, that quant_below_idx is nondecreasing, and
 * that gamma and quant_below_idx have the same length.
 *
 * -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
 *
 * The quantile calculation part of this function is an implementation of the R
 * function stats::quantile under the default choice of empirical quantile.  See
 * the help page for stats::quantile for more details (the default there is type
 * 7).  The formula is (for n the number of observations in the data):
 *
 *     F^{-1}(prob) = ((1 - gamma) * x_j) + (gamma * x_{j+1})
 *
 * where
 *
 *     j     = floor(prob * (n - 1))
 *     gamma = (prob * (n - 1)) - floor(prob * (n - 1))
 */

arma::vec quantile(double* values,
		   int data_len,
		   arma::Mat<int>& quant_below_idx, 
		   arma::vec& gamma) {

    int nprobs = quant_below_idx.n_elem;  // num empirical quantiles
    int curr_size = data_len;            // num samples in data
    int relative_below_idx[nprobs];       // index relative to previous index

    double* curr_start;  // start of the current subset of the data
    double* end;         // 1 past the end of the data
    double* kplus1;      // 1 past the start of the current subset of the data
    double* curr;        // pointer to the current value
    double* curr_min;    // pointer to the current minimum

    int below_idx;  // index leq to desired index location
    int above_idx;  // index one greater than below_idx(k)
    double temp;    // used for data swap

    arma::vec out(nprobs);  // storage for empirical quantiles

    // Calculate initial starting location and absolute end of data
    curr_start = values;
    end = curr_start + data_len;

    // Calculate indices relative to their previous index
    relative_below_idx[0] = quant_below_idx(0);
    for (int i = 1; i < nprobs; i++) {
	relative_below_idx[i] = quant_below_idx(i) - quant_below_idx(i - 1);
    }

    /* Each iteration calculates the empirical quantile value of the k-th
     * probability.  The main task is to partially sort the data so that the
     * values in the quant_below_idx(k)-th index and the index immediately
     * following are the same is if the data had been fully sorted.
     *
     * This is done by partially sorting the data with indices higher than the
     * partial sort in the previous iteration.
     */ 
    for (int k = 0; k < nprobs; k++) {
	    
	/* Partial sorting on k-th quantile index; i.e. data is permuted so that
	 * every value with index less than k-th quantile index of the data is
	 * smaller than value at the k-th quantile index, and every value with
	 * index larger than k-th quantile index of the data is larger.
	 */
	rPsort(curr_start, curr_size, relative_below_idx[k]);

	/* Find the index of the (k + 2)-th largest value.  We start
	 * curr_min_ptr at the (k + 1)-th index b/c we know every value with
	 * smaller index is <= to the value there.
	 */ 
	kplus1 = curr_start + relative_below_idx[k] + 1;
	curr_min = kplus1;
	for (curr = curr_min; curr < end; curr++) {
	    if (*curr < *curr_min) {
		curr_min = curr;
	    }
	}

	// Swap the value in the (k + 1)-th index with the (k + 2)-th largest value
	temp = *kplus1;
	*kplus1 = *curr_min;
	*curr_min = temp;
		
	/* Calculate empirical quantile value.  Note that we transpose the
	 * indices so as to save the quantile data in variable by quantile
	 * level shape.  It is (k + 1) rather than k b/c we fill the 0-th
	 * column with the true values of beta.
	 */
	below_idx = quant_below_idx(k);
	above_idx = below_idx + 1;
	out(k) = ((1 - gamma(k)) * values[below_idx]) + (gamma(k) * values[above_idx]);

	// Reflect the subset of the data that needs partial sorting
	curr_start += relative_below_idx[k];
	curr_size = end - curr_start;
    }

    return out;
}




/* Calculate the empirical quantiles for each column of X and quantile levels
 * given by prob, and return as elements of an arma::mat.  The first column of
 * the return object is the true values of the variables, and the (j, k + 1)
 * element of the return object is the quantile value of the j-th variable at
 * the k-th quantile level.
 *
 * PRE: assumes that true_beta is the same length as the number of columns of X,
 * that X has >= 2 rows of data, and that prob is nondecreasing and has values
 * in (0, 1).
 */

arma::mat quantile_table(arma::vec& true_beta, arma::mat& X, arma::vec& prob) {

    int p = X.n_cols;          // number of variables
    int nprobs = prob.n_elem;  // number of probabilities
    int nsamp = X.n_rows;      // number of samples
    double* col_start;         // first element in X.col(j)

    // Container for quantile values
    arma::mat out(nprobs + 1, p);


    // Quantile indices, interpolation vals ----------------

    double desired;  // idx you might choose point if you had inf points
    arma::Col<int> below_idx(nprobs);  // closest index below desired point
    arma::vec gamma(nprobs);           // linear interpolation value

    /* Calculate data indices at which the quantiles are calculated and the
     * corresponding linear interpolation values
     */
    for (int k = 0; k < nprobs; k++) {
	// Minus 1 b/c of 0-based indexing
	desired = prob(k) * (nsamp - 1);
	below_idx(k) = (int) desired;
	gamma(k) = desired - below_idx(k);
    }


    // Fill in values of table -----------------------------

    // Fill in true values of beta in first row
    out.row(0) = true_beta.t();

    /* Each iteration fills in rows 1-3 of column j with the empirical quantile
     * levels of the j-th column of X
     */
    col_start = X.memptr();
    for (int j = 0; j < p; j++) {

	// Fill in empirical quantile levels
	out.submat(1, j, 3, j) = quantile(col_start, nsamp, below_idx, gamma);

	// Update col_start to point to first element of (j + 1)-th column
	col_start += nsamp;
    }

    return out.t();
}
