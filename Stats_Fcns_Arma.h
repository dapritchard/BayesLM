
arma::vec mvrnorm_chol(arma::vec& m, arma::mat& V);
arma::vec mvrnorm_eigen(arma::vec& m, arma::mat& V);
arma::vec quantile(arma::mat& X, double prob);
arma::vec sample_beta(int p, double prop_nonzero, double true_beta_sd);

