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


#n <- 50
#n <- 1000
n <- 100
#n <- 200
q <- 0.2
#q <- 0.1
alpha <- c(
    0.000001,
    0.000005,
    0.00001,
    0.00005,
    0.0001,
    0.0005,
    0.001,
    0.005,
    0.01,
    0.05)
alpha <- c(alpha, 1:9/10)

pmfs <- allpmfs(n, q)
pr <- privratio(pmfs)
quants <- pmfquantile(pmfs, alpha, interpolate = TRUE)


## For a given collection size n, plot the privacy ratio as a function of the
## number of 1s in the synthetic collection s = 0:n, for each number of 1s in
## the original collection m = 1:n.
## Optionally indicate the values at a specific quantiles of A(m) for each m.
## Privacy ratio and quantile values should be supplied as tables of
## (m, s, pr) and (m, alpha, qs) respectively.
## If 'segments' is TRUE, connect dots across m values for each alpha.
## If 'maxovers' is TRUE, plot only the largest pr value if there are mutliple
## for a given s value and alpha.
plotpr <- function(pr, quants = NULL, segments = FALSE, maxovers = FALSE) {
    n <- pr[, max(s)]
    ## Create the basic plot with separate lines for each m.
    gp <- ggplot(pr, aes(x = s, y = pr, group = factor(m))) +
        geom_line(colour = "grey") +
        ylim(pr[, range(pr)]) +
        theme(text = element_text(size = 10))
    if(!is.null(alpha)) {
        setkey(pr, m, s)
        setkey(quants, m, qs)
        ## Look up the privacy ratio values at each s point.
        dd <- pr[quants, nomatch = 0]
        ## Optionally only plot the max if there are multiple y values for
        ## given s.
        if(maxovers) dd <- dd[, list(pr = max(pr)), by = list(s, alpha)]
        gp <- gp + geom_point(data = dd, aes(colour = factor(alpha)))
        ## Add line segments to connect the points if required.
        if(segments) {
            ddseg <- dd[, list(x = s[-.N], xend = s[-1], y = pr[-.N],
                yend = pr[-1]), by = alpha]
            gp <- gp + geom_segment(data = ddseg, aes(x = x, xend = xend,
                y = y, yend = yend, colour = factor(alpha)))
        }
        ## Add indicators for the max privacy ratio value for each alpha.
        ddmax <- dd[, .SD[which.max(pr)], by = alpha]
        gp <- gp + geom_point(data = ddmax) +
            geom_point(data = ddmax, aes(colour = factor(alpha)), size = 0.5)
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

pdf("pr-quantiles.pdf", 12, 8)
plotpr(pr, quants, TRUE, TRUE)
dev.off()


## Make plots for multiple values of n and q.
#nvals <- c(100, 500)
#qvals <- c(0.1, 0.2, 0.4)
#alpha <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9,
#    0.00001, 0.00005, 0.000001, 0.000005)
#qnplots <- setNames(lapply(qvals, function(q) {
#    setNames(lapply(nvals, function(n) { plotpr(n, q, alpha) }),
#        as.character(nvals))
#}), as.character(qvals))
#
#pdf("ratio-at-quantiles.pdf", 16, 8)
#lapply(qnplots, function(pp) { multiplot(plotlist = pp, cols = 2) })
#dev.off()
#

## Make plots of pmfs indicating quantiles.
pmfA <- allpmfs(50, 0.2)
dpmf <- rbindlist(mapply(function(i, v) {
    data.table(m = i - 1, s = seq_along(v) - 1, pmf = v)
}, seq_along(pmfA), pmfA, SIMPLIFY = FALSE, USE.NAMES = FALSE))
dpmf[, cdf := cumsum(pmf), by = m]

qA <- pmfquantile(pmfA, alpha = c(0.01, 0.05, 0.1))
setkey(qA, m, qs)
setkey(dpmf, m, s)
qA <- dpmf[qA]

#pvar <- "pmf"
pvar <- "cdf" 
pvarn <- as.name(pvar)

gp <- ggplot(dpmf, aes(x = s, y = eval(pvarn))) + geom_line() +
    geom_point(size = 0.5) +
    facet_wrap(~ m) +
    theme(text = element_text(size = 6)) +
    labs(title = sprintf("%s of A(m) over m = 0:n with quantiles indicated",
        toupper(pvar)),
        x = "x", y = "Prob", colour = "Quantile prob")
gp <- gp + geom_point(data = qA, aes(colour = factor(alpha))) +
    geom_hline(data = qA, aes(yintercept = eval(pvarn), colour = factor(alpha)))
if(pvar == "cdf") gp <- gp + ylim(0, 0.2)

pdf(sprintf("%ss-quantiles.pdf", pvar), 12, 10); gp; dev.off()


