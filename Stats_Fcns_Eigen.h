
Eigen::VectorXd mvrnorm_chol(Eigen::VectorXd& m, Eigen::MatrixXd& V);
Eigen::VectorXd mvrnorm_eigen(Eigen::VectorXd& m, Eigen::MatrixXd& V);

Eigen::MatrixXd quantile_table(Eigen::VectorXd& true_beta, 
			       Eigen::MatrixXd& X, 
			       Eigen::VectorXd& prob);

Eigen::MatrixXd quantile_table(double true_beta, 
			       Eigen::VectorXd& X, 
			       Eigen::VectorXd& prob);

Eigen::VectorXd quantile(double* values,
			 int data_len,
			 Eigen::VectorXi& quant_below_idx, 
			 Eigen::VectorXd& gamma);

Eigen::VectorXd sample_beta(int p, double prop_nonzero, double true_beta_sd);
Eigen::MatrixXd matr_randn(int nrow, int ncol, double sd);
