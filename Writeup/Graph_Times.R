
library(calibrate)
library(magrittr)
expression(paste("Plot of ", alpha^beta, " versus ", hat(mu)[0])))
xvals <- 1:3
xlab <- "Data size"
ylab <- expression(paste(log[2], " time in secs"))
titlemain <- "Graph as p increases"
titlevec <- c("inv time", "mvn time", "sampler time", "total time")
colorvec <- c("red", "blue", "green", "purple", "orange")





# Expands graph boundaries by (100 * prop_expand)%
expand_lim <- function(lim_vals, prop_expand) {
    expand_by <- (lim_vals[2L] - lim_vals[1L]) * prop_expand / 2
    c(lim_vals[1L] - expand_by, lim_vals[2L] + expand_by)
}



graph_time_type <- function(time_dat, xvals, type_nm, xlab_nm, ylab_nm, title_nm, colorvec) {

    # Restrict data to the time type of interest
    time_type_dat <- subset(time_dat, time_type == type_nm)
    program_nm <- unique(time_type_dat$program)

    # Calculate graph boundaries
    allvals <- time_type_dat$time
    xlim <- c(min(xvals), max(xvals)) %>% expand_lim(., 0.25)
    ylim <- c(min(allvals), max(allvals)) %>% expand_lim(., 0.1)

    # Plot without any points
    plot(-100, xlim=xlim, ylim=ylim, xlab=xlab_nm, ylab=ylab_nm, main=title_nm)

    # Each iteration draws the lines for one program
    for (i in seq_along(program_nm)) {

        # Plot the lines for program_nm[i]
        yvals <- time_type_dat$time[ time_type_dat$program == program_nm[i] ]
        lines(xvals, yvals, type="b", pch=20, col=colorvec[i])

        # Add points to the lines
        y_fmt <- format(round(yvals, 2), nsmall=2)
        calibrate::textxy(X      = xvals,
                          Y      = yvals,
                          lab    = y_fmt,
                          m      = c(mean(xlim), mean(ylim)),
                          cex    = 1,
                          offset = 0.75,
                          col    = colorvec[i])
    }
}



graph_times <- function(time_obs, xvals, xlab_nm, ylab_nm, titlevec, titlemain, colorvec) {

    # Set up the layout of the graph
    par(oma=c(0, 0, 3, 0), las=1, family="serif")
    layout(matrix(1:9, nrow=3, ncol=3, byrow=TRUE),
           widths = c(0.475, .015, .475),
           heights = c(0.475, .05, .475))

    # Margins for each plot
    margin_vals <- c(4.1, 4, 1.5, 0.1)
    zerovec <- rep(0, 4)

    # 1, 1
    par(mar=margin_vals)
    graph_time_type(time_obs, xvals, "inverse", xlab_nm, ylab_nm, titlevec[1L], colorvec)

    # 1, 2
    par(mar=zerovec)
    frame()

    # 1, 3
    par(mar=margin_vals)
    graph_time_type(time_obs, xvals, "mvnsamp", xlab_nm, ylab_nm, titlevec[2L], colorvec)

    # 2, 1
    par(mar=zerovec)
    frame()

    # 2, 2
    plot(-100, -100, bty="n", xlim=c(0, 1), ylim=c(0, 1), axes=FALSE, xlab=NULL)
    legend(x=-1.35, y=2.5, legend=c("R only", "Armadillo", "Eigen"),  lty=c(1, 1, 1),
           lwd=c(1, 1, 1), col=colorvec[c(1, 2, 3)], xpd=NA, box.lwd=0.5)

    # 2, 3
    par(mar=zerovec)
    frame()

    # 3, 1
    par(mar=margin_vals)
    graph_time_type(time_obs, xvals, "mat_ops", xlab_nm, ylab_nm, titlevec[3L], colorvec)

    # 3, 2
    par(mar=zerovec)
    frame()

    # 3, 3
    par(mar=margin_vals)
    graph_time_type(time_obs, xvals, "total", xlab_nm, ylab_nm, titlevec[4L], colorvec)

    mtext(titlemain, outer=TRUE, cex = 1.25, padj=-1)
}
