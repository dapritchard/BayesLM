
# Collects time observations for various functions and levels of n, p, nsamp,
# and decomposition method.  fcn_nm_vec is a character vector of function names
# corresponding to functions to call. nvec and pvec are numeric vectors of equal
# length providing the levels of n and p to run each function with.  nrepl,
# nsamp and decmp are single values that stay constant for every function call.
#
# Each function is run with every n, p level as given by the k-th element of
# nvec and pvec for every k in 1, ..., K.  In other words there are K levels of
# n, p, and not K^2 levels.

collect_time_data <- function(fcn_nm_vec, nvec, pvec, nrepl, nsamp=1e4, decmp="chol") {

    out <- list()
    ctr <- 1L

    # Each iteration sets of of the functions named in fcn_nm_vec and collects
    # data for each of the (n, p, nsamp) triplets, appending one element to out
    # for each triplet
    for (fcn_nm in fcn_nm_vec) {
        fcn <- get(fcn_nm)

        # Each iteration collects data for a single function named in fcn_nm_vec
        # and a single set of (n, p, nsamp) values and appends an element to out
        for (i in seq_along(nvec)) {
            n <- nvec[i]
            p <- pvec[i]
            datavals <- replicate(nrepl, fcn(n=n, p=p, nsamp=nsamp, decomp_method=decmp))

            out[[ctr]] <- list(fcn_nm = fcn_nm,
                               n      = n,
                               p      = p,
                               dat    = datavals)
            ctr <- ctr + 1L
        }
    }

    list(nvec  = nvec,
         pvec  = pvec,
         nrepl = nrepl,
         nsamp = nsamp,
         decmp = decmp,
         time  = out)
}




# Formats the data provided by collect_time_data into a list containing a
# data.frame of the mean observed runtimes

format_time_data <- function(time_data) {

    datlen <- length(time_data$time)
    tablen <- 4 * datlen

    program   <- vector("character", tablen)
    time_type <- rep(c("inverse", "mvnsamp", "total", "mat_ops"), datlen)
    n         <- vector("numeric", tablen)
    p         <- vector("numeric", tablen)
    time_vals <- vector("numeric", tablen)

    ctr <- 1L

    # Each iteration adds 4 times values and the corresponding argument
    # specifications used in the function call to the data vectors
    for (elem in time_data$time) {

        idx <- ctr:(ctr + 3L)
        mean_time <- apply(elem$dat, 1, mean)

        program[idx] <- elem$fcn_nm
        n[idx] <- elem$n
        p[idx] <- elem$p

        time_vals[ idx[-4L] ] <- mean_time
        time_vals[ idx[4L] ] <- mean_time[3L] - sum( mean_time[1:2] )

        ctr <- ctr + 4L
    }

    # Construct time values data.frame
    out_time <- data.frame(program   = program,
                           time_type = time_type,
                           n         = n,
                           p         = p,
                           time      = log2(time_vals),
                           # options
                           stringsAsFactors=FALSE)

    # Sort data.frame by program and time_type
    time_type <- factor(time_type, c("inverse", "mvnsamp", "mat_ops", "total"))
    sort_time <- out_time[order(program, time_type), ]

    # Return results
    list(nvec  = time_data$nvec,
         pvec  = time_data$pvec,
         nrepl = time_data$nrepl,
         nsamp = time_data$nsamp,
         decmp = time_data$decmp,
         dat   = sort_time)
}
