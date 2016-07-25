
Eigen::VectorXd mvrnorm_chol(Eigen::VectorXd& m, Eigen::MatrixXd& V);
Eigen::VectorXd mvrnorm_eigen(Eigen::VectorXd& m, Eigen::MatrixXd& V);
Eigen::VectorXd quantile(Eigen::MatrixXd& X, double prob);
double quantile(Eigen::VectorXd& X, double prob);
Eigen::VectorXd sample_beta(int p, double prop_nonzero, double true_beta_sd);
Eigen::MatrixXd matr_randn(int nrow, int ncol, double sd);
