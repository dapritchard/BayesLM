
source("Bayes_LM.R")
source("Bayes_LM_Rcpp_Arma.R")
source("Bayes_LM_Rcpp_Eigen.R")
source("Writeup/Collect_Time_Data.R")
source("Writeup/Graph_Times.R")

fcns_cpp_nm <- c("bayes_lm_rcpp_arma", "bayes_lm_rcpp_eigen")
fcns_all_nm <- c("bayes_lm_r", "bayes_lm_rcpp_arma", "bayes_lm_rcpp_eigen")
fcns_open_nm <- c("bayes_lm_r", "bayes_lm_rcpp_arma")


# As n --> infty for R + Atlas, Arma + Atlas, and Eigen
if (! file.exists("Writeup/data/n_increase_atlas.RData")) {

    n_increase_atlas_cpp <- collect_time_data(fcn_nm_vec = fcns_cpp_nm,
                                               nvec      = 2^c(6, 9, 12, 15, 18),
                                               pvec      = rep(100, 5),
                                               nrepl     = 10,
                                               nsamp     = 1e4)

    n_increase_atlas_r <- collect_time_data(fcn_nm_vec = "bayes_lm_r",
                                            nvec       = 2^c(6, 9, 12, 15, 18),
                                            pvec       = rep(100, 5),
                                            nrepl      = 5,
                                            nsamp      = 1e4)

    save(n_increase_atlas_cpp, n_increase_atlas_r, file="Writeup/data/n_increase_atlas.RData")
}


# As n --> infty for R + OpenBLAS and Arma + OpenBLAS
if (! file.exists("Writeup/data/n_increase_openBLAS.RData")) {

    n_increase_openBLAS <- collect_time_data(fcn_nm_vec = fcns_open_nm,
                                             nvec       = 2^c(6, 9, 12, 15, 18),
                                             pvec       = rep(100, 5),
                                             nrepl      = 10,
                                             nsamp      = 1e4)

    save(n_increase_openBLAS, file="Writeup/data/n_increase_openBLAS.RData")
}


# As p --> infty for R + OpenBLAS
if (! file.exists("Writeup/data/p_increase_r_openblas.RData")) {

    p_increase_r_openblas <- collect_time_data(fcn_nm_vec = "bayes_lm_r",
                                               nvec       = rep(2000, 5),
                                               pvec       = 2^(6:10),
                                               nrepl      = 6,
                                               nsamp      = 1e4)

    save(p_increase_r_openblas, file="Writeup/data/p_increase_r_openblas.RData")
}


# As p --> infty for Arma + OpenBLAS and Eigen
if (! file.exists("Writeup/data/p_increase_arma_openblas_eigen.RData")) {

    p_increase_arma_openblas_eigen <- collect_time_data(fcn_nm_vec = fcns_cpp_nm,
                                                        nvec       = rep(2000, 5),
                                                        pvec       = 2^(6:10),
                                                        nrepl      = 10,
                                                        nsamp      = 1e4)

    save(p_increase_arma_openblas_eigen, file="Writeup/data/p_increase_arma_openblas_eigen.RData")
}


# As p --> infty for R + ATLAS
if (! file.exists("Writeup/data/p_increase_r_atlas.RData")) {

    p_increase_r_atlas <- collect_time_data(fcn_nm_vec = "bayes_lm_r",
                                            nvec       = rep(2000, 5),
                                            pvec       = 2^(6:10),
                                            nrepl      = 6,
                                            nsamp      = 1e4)

    save(p_increase_r_atlas, file="Writeup/data/p_increase_r_atlas.RData")
}


# As p --> infty for Arma + ATLAS
if (! file.exists("Writeup/data/p_increase_arma_atlas.RData")) {

    p_increase_arma_atlas <- collect_time_data(fcn_nm_vec = "bayes_lm_rcpp_arma",
                                               nvec       = rep(2000, 5),
                                               pvec       = 2^(6:10),
                                               nrepl      = 6,
                                               nsamp      = 1e4)

    save(p_increase_arma_atlas, file="Writeup/data/p_increase_arma_atlas.RData")
}



# x <- format_time_data(p_increase_r_atlas)$dat





# Graph n --> infty  -----------------------------------------------------------

# Pull in data
load("Writeup/data/n_increase_atlas.RData")
load("Writeup/data/n_increase_openBLAS.RData")

# Mung data into data.frame form (after taking averages)
obs_cpp_atlas <- format_time_data(n_increase_atlas_cpp)
obs_r_atlas <- format_time_data(n_increase_atlas_r)
obs_openBLAS <- format_time_data(n_increase_openBLAS)

# Rename the program name for OpenBLAS so that we can distinguish from Atlas obs
obs_openBLAS$dat$program <- paste0(obs_openBLAS$dat$program, "_openBLAS")

# Combine into a single data.frame
n_increase <- rbind(obs_cpp_atlas$dat, obs_r_atlas$dat, obs_openBLAS$dat, stringsAsFactors=FALSE)

# Order in which functions will be graphed
sort(unique(n_increase$program))

# Construct graph parameter specifications
xvals <- c(6, 9, 12, 15, 18)
xlab <- expression(paste(log[2], "  number of observations"))
ylab <- expression(paste(log[2], "  time in secs"))
titlevec <- c("Inverse", "Sampling MVN", "Matrix Operations", "Total")
titlemain <- "Graph as n increases (p=100, samples=10,000)"
colorvec <- c(r_atlas="red", arma_atlas="blue", arma_openblas="purple",
              eigen="green2", r_openblas="orange")
legendvec <- c("R + ATLAS", "Arma + ATLAS", "Arma + OpenBLAS", "Eigen", "R + OpenBLAS")

# Build graph
pdf("Writeup/figure/n_increase.pdf", width=8, height=6)
graph_times(n_increase, xvals, xlab, ylab, titlevec, titlemain, colorvec, legendvec)
dev.off()





# Graph p --> infty  -----------------------------------------------------------

# Pull in data
load("Writeup/data/p_increase_arma_atlas.RData")
load("Writeup/data/p_increase_arma_openblas_eigen.RData")
load("Writeup/data/p_increase_r_atlas.RData")
load("Writeup/data/p_increase_r_openblas.RData")

# Mung data into data.frame form (after taking averages)
obs_arma_atlas <- format_time_data(p_increase_arma_atlas)
obs_arma_openblas_eigen <- format_time_data(p_increase_arma_openblas_eigen)
obs_r_atlas <- format_time_data(p_increase_r_atlas)
obs_r_openblas <- format_time_data(p_increase_r_openblas)

# Rename the program name so that we can distinguish from each other
obs_arma_atlas$dat$program <- paste0(obs_arma_atlas$dat$program, "_atlas")
obs_arma_openblas_eigen$dat$program[
    obs_arma_openblas_eigen$dat$program == "bayes_lm_rcpp_arma"
] <- "bayes_lm_rcpp_arma_openblas"
obs_r_atlas$dat$program <- paste0(obs_r_atlas$dat$program, "_atlas")
obs_r_openblas$dat$program <- paste0(obs_r_openblas$dat$program, "_openblas")

# Combine into a single data.frame
p_increase <- rbind(obs_arma_atlas$dat,
                    obs_arma_openblas_eigen$dat,
                    obs_r_atlas$dat,
                    obs_r_openblas$dat,
                    # options
                    stringsAsFactors=FALSE)

# Order in which functions will be graphed
sort(unique(p_increase$program))

# Construct graph parameter specifications
xvals <- 6:10
xlab <- expression(paste(log[2], "  number of variables"))
ylab <- expression(paste(log[2], "  time in secs"))
titlevec <- c("Inverse", "Sampling MVN", "Matrix Operations", "Total")
titlemain <- "Graph as p increases (n=2,000, samples=10,000)"
colorvec <- c(r_atlas="red", arma_atlas="blue", arma_openblas="purple",
              eigen="green2", r_openblas="orange")
legendvec <- c("R + ATLAS", "Arma + ATLAS", "Arma + OpenBLAS", "Eigen", "R + OpenBLAS")

# Build graph
pdf("Writeup/figure/p_increase.pdf", width=8, height=6)
graph_times(p_increase, xvals, xlab, ylab, titlevec, titlemain, colorvec, legendvec)
dev.off()
