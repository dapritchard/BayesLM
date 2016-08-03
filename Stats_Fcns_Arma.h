
arma::vec mvrnorm_chol(arma::vec& m, arma::mat& V);
arma::vec mvrnorm_eigen(arma::vec& m, arma::mat& V);

arma::mat quantile_table(arma::vec& true_beta, arma::mat& X, arma::vec& prob);
arma::vec quantile(double* values,
		   int data_len,
		   arma::Mat<int>& quant_below_idx, 
		   arma::vec& gamma);


/* arma::vec quantile(arma::mat& X, double prob); */
arma::vec sample_beta(int p, double prop_nonzero, double true_beta_sd);

