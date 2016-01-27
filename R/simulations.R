#############################################################################
###
###  Numerical simulations and tests for k-randomization project.
###
#############################################################################


library(data.table)
library(ggplot2)

## multiplot() function used below can be found at
## https://github.com/dzeber/work-tools/blob/master/R/utils/ggplot.R.

source("k-rand-lib.R")

## For a given collection size n, plot the privacy ratio as a function of the
## number of 1s in the synthetic collection s = 0:n, for each number of 1s in
## the original collection m = 1:n.
## Optionally indicate the values at a specific quantiles of A(m) for each m.
plotpr <- function(n, q = 0.2, alpha = NULL) {
    ## First compute the table of privacy ratio values.
    pr <- privratio(n, q)
    pr <- rbindlist(mapply(function(m, v) {
        data.table(m = m, s = seq_along(v) - 1, pr = v)
    }, seq_along(pr), pr, USE.NAMES = FALSE, SIMPLIFY = FALSE))
    ## Create the basic plot with separate lines for each m.
    gp <- ggplot(pr, aes(x = s, y = pr, group = factor(m))) +
        geom_line(colour = "grey") +
        ylim(q/(1-q), (1-q)/q) +
        theme(text = element_text(size = 10))
    if(!is.null(alpha)) {
        setkey(pr, m, s)
        ## Compute quantiles for each probability value alpha.
        qvals <- rbindlist(lapply(alpha, function(a) {
            v <- as.numeric(lapply(1:n, qbinomsum, n, q, a))
            data.table(m = 1:n, s = v, alpha = a)
        }))
        ## Look up the privacy ratio values at each s point.
        qvals <- pr[qvals]
        gp <- gp + geom_point(data = qvals, aes(colour = factor(alpha)))
        ## Add indicators for the max privacy ratio value for each alpha.
        maxqvals <- qvals[, .SD[which.max(pr)], by = alpha]
        gp <- gp + geom_point(data = maxqvals)
    }
    gp <- gp + labs(
        title = paste("Privacy ratio vs. s plotted separately for each m",
            sprintf("(n = %s, q = %s)", n, q),
            if(!is.null(alpha)) paste(
                "\nwith indicators at quantile values for each m",
                "(max indicated in black)")
            else ""))
    if(!is.null(alpha)) gp <- gp + labs(colour = "Quantile prob")
    gp
}


## Make plots for multiple values of n and q.
nvals <- c(100, 500)
qvals <- c(0.1, 0.2, 0.4)
alpha <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9,
    0.00001, 0.00005, 0.000001, 0.000005)
qnplots <- setNames(lapply(qvals, function(q) {
    setNames(lapply(nvals, function(n) { plotpr(n, q, alpha) }),
        as.character(nvals))
}), as.character(qvals))

pdf("ratio-at-quantiles.pdf", 16, 8)
lapply(qnplots, function(pp) { multiplot(plotlist = pp, cols = 2) })
dev.off()

