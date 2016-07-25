
library(calibrate)
library(magrittr)

xvals <- 1:5
xlab <- "Data size"
ylab <- "log2 time in secs"
titlemain <- "Graph as p increases"
titlevec <- c("inv time", "mvn time", "sampler time", "total time")
colorvec <- c("red", "blue", "green", "purple", "orange")

tinv <- lapply(1:5, function(x) rnorm(5) + x)
tmvn <- tsam <- ttot <- tinv


graph_times <- function(ctime, xvals, xlab, ylab, titlevec, titlemain, colorvec) {

    time_graph <- function(x, ylist, xlab, ylab, title_nm, colorvec) {

        expand_lim <- function(lim_vals, prop_expand) {
            expand_by <- (lim_vals[2] - lim_vals[1]) * prop_expand / 2
            c(lim_vals[1] - expand_by, lim_vals[2] + expand_by)
        }

        allvals <- unlist(ylist)
        xlim <- c(min(x), max(x)) %>% expand_lim(., 0.25)
        ylim <- c(min(allvals), max(allvals)) %>% expand_lim(., 0.1)

        plot(-100, xlim=xlim, main=title_nm, ylim=ylim, xlab=xlab, ylab=ylab)
        for (i in 1:length(ylist)) {
            lines(x, ylist[[i]], type="b", pch=20, col=colorvec[i])
            sani_labels <- format(round(ylist[[i]], 2), nsmall=2)
            calibrate::textxy(X      = x,
                              Y      = ylist[[i]],
                              lab    = sani_labels,
                              m      = c(mean(xlim), mean(ylim)),
                              cex    = 1,
                              offset = 0.75,
                              col    = colorvec[i])
        }
    }

    par(oma=c(0, 0, 3, 0), las=1, family="serif")
    par(las=1)
    layout(matrix(1:9, nrow=3, ncol=3, byrow=TRUE),
          widths = c(0.475, .015, .475),
          heights = c(0.475, .05, .475))

    # 1, 1
    par(mar=c(4, 4, 1.5, 0.1))
    time_graph(xvals, ctime$inv, xlab, ylab, titlevec[1], colorvec)

    # 2, 1
    par(mar=rep(0, 4))
    frame()

    # 3, 1
    par(mar=c(4, 4, 1.5, 0.1))
    time_graph(xvals, ctime$mvn, xlab, ylab, titlevec[2], colorvec)

    # 1, 2 / 2, 2 / 3, 2
    par(mar=rep(0, 4))
    frame()
    frame()
    frame()

    # 1, 3
    par(mar=c(4, 4, 1.5, 0.1))
    time_graph(xvals, ctime$sam, xlab, ylab, titlevec[3], colorvec)

    # 2, 3
    par(mar=rep(0, 4))
    frame()

    # 3, 3
    par(mar=c(4, 4, 1.5, 0.1))
    time_graph(xvals, ctime$tot, xlab, ylab, titlevec[4], colorvec)

    mtext(titlemain, outer=TRUE, cex = 1.25, padj=-1)
}




transform_times <- function(time_list) {
    library(magrittr)

    # Set up list structure.  The first level is the time being measured.
    retlist <- setNames(vector("list", 4), c("inv", "mvn", "tot", "sam"))

    # Second level of list is the function / program being used
    for (k in 1:4) {
        retlist[[ k ]] <- setNames(vector("list", 5),
                               c("ronly", "arma", "eigen", "rcpp_arma", "rcpp_eigen"))

        # Third level of list is the elapsed times as the data changes
        for (i in 1:5) {
            retlist[[ k ]][[ i ]] <- vector("numeric", 5)
        }
    }

    # i steps through the function / program being used
    for (i in 1:5) {

        # j steps through the size of the data being used
        for (j in 1:5) {

            # k steps through the type of time being measured
            for (k in 1:3) {
                retlist[[ k ]][[ i ]][ j ] <- time_list[[ i ]][[ j ]][ k, ] %>% mean %>% log2
            }
            retlist[[ 4 ]][[ i ]][[ j ]] <- ( time_list[[ i ]][[ j ]][ 3, ] -
                                          time_list[[ i ]][[ j ]][ 1, ] -
                                          time_list[[ i ]][[ j ]][ 2, ] ) %>% mean %>% log2
        }
    }

    retlist
}
