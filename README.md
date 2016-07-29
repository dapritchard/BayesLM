
## A case study in the implementation of a Bayesian linear model using R and / or C++ libraries

As an exercise I decided to implement a Bayesian linear model using the `C++`
libraries `Armadillo` and `Eigen`, and then to interface the `C++` code from `R`
via `Rcpp`. The main purpose of this exercise was to obtain some familiarity
with these libraries and packages.  Additionally I am also interested in
comparing the speeds of the implementations for this problem and practicing
profiling techniques on the software.  

A more detailed explanation of the functions / programs found in this repository
can be found at `Writeup/Bayes_LM_Writeup.pdf`. The file also displays my
findings and experiences in terms of speed comparisons and program profiling.  

The root directory contains the following files:

1. `bayes_lm_r`: an `R` function
    * `Bayes_LM.R`
    * `Check_Valid_Input.R`
	
2. `bayes_lm_arma`: a `C++`-only implementation using the `Armadillo` library (an executable)
    * `Bayes_LM_Arma.cpp`
    * `Parse_Args.cpp`
    * `Stats_Fcns_Arma.cpp`
  
3. `bayes_lm_eigen`: a `C++`-only implementation using the `Eigen` library (an executable)
    * `Bayes_LM_Eigen.R`
    * `Parse_Args.cpp`
    * `Stats_Fcns_Eigen.cpp`
  
4. `bayes_lm_rcpp_arma`: an `R` function internally calling a workhorse `C++` function constructed using the `Armadillo` library
    * `Bayes_LM_Rcpp_Arma.R`
    * `Check_Valid_Input.R`
    * `Stats_Fcns_Arma.cpp`
  
5. `bayes_lm_rcpp_eigen`: an R function internally calling a workhorse `C++` function constructed using the `Eigen` library
    * `Bayes_LM_Rcpp_Eigen.R`
    * `Check_Valid_Input.R`
    * `Stats_Fcns_Eigen.cpp`

![n increases](Writeup/figure/n_increase.pdf "Title")
