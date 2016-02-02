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

n <- 50
#n <- 1000
#n <- 100
#n <- 200
q <- 0.2
#q <- 0.1
alpha <- c(
#    0.000001,
#    0.000005,
#    0.00001,
#    0.00005,
#    0.0001,
#    0.0005,
#    0.001,
#    0.005,
    0.01,
    0.05,
    0.1)
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

pdf("pr-quantiles-n50.pdf", 12, 8)
plotpr(pr, quants, TRUE, FALSE)
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


##########################################################
##########################################################

## Make plots to explore distributions and quantiles.
#n <- 50
n <- 100
q <- 0.2
#alpha <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.9)
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
    0.05,
    0.1,
    0.99,
    0.999,
    0.9999
    )
alpha <- c(alpha, 1:9/10)

pmfs <- allpmfs(n, q)
quants <- pmfquantile(pmfs, alpha, interpolate = TRUE)
pr <- privratio(pmfs)
## Privacy ratio at quantiles.
setkey(pr, m, s)
setkey(quants, m, qs)
prquants <- quants[pr][, list(m, alpha, qs, pr, qs.interp, interp)]
prquants[, pr.interval := rev(diff(c(rev(pr), 0))), by = m]
prquants[, pr.interp := pr + interp * pr.interval][, pr.interval := NULL]
prquants <- prquants[!is.na(alpha)]
## Prob ratios at quantiles.
setkey(ddist, m, s)
ddistquants <- quants[ddist][, list(m, alpha, qs, qs.interp, interp,
    pmf = i.pmf, cdf = i.cdf)]
ddistquants[, probratio := 

## Create table of distribution.
ddist <- rbindlist(mapply(function(i, v) {
    data.table(m = i - 1, s = seq_along(v) - 1, pmf = v)
}, seq_along(pmfs), pmfs, SIMPLIFY = FALSE, USE.NAMES = FALSE))
ddist[, cdf := cumsum(pmf), by = m]

## Plot pmfs by m with quantiles indicated.
plotpmfs <- function(ddist, quants = NULL) {
    gp <- ggplot(ddist, aes(x = s, y = pmf)) + geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ m) +
        theme(text = element_text(size = 6)) +
        labs(title = sprintf("PMF of A(m) over m = 0:%s", ddist[, max(m)]),
            x = "s", y = "Prob")
    if(!is.null(quants)) {
        gp <- gp + geom_point(data = quants,
                aes(x = qs, colour = factor(alpha))) +
            ## Indicate value at quantile with horizontal line.
            geom_hline(data = quants,
                aes(yintercept = pmf, colour = factor(alpha))) +
            labs(title = sprintf("%s with quantiles indicated", gp$labels$title),
               colour = "Quantile prob")
    }
    gp
}

## Plot cdfs by m with quantiles indicated.
plotcdfs <- function(ddist, quants = NULL) {
    gp <- ggplot(ddist, aes(x = s, y = cdf)) + geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ m) +
        theme(text = element_text(size = 6)) +
        labs(title = sprintf("CDF of A(m) over m = 0:%s", ddist[, max(m)]),
            x = "s", y = "Prob")
    if(!is.null(quants)) {
        gp <- gp + geom_point(data = quants,
                aes(x = qs, colour = factor(alpha))) +
            ## Indicate value at quantile with horizontal line.
            geom_hline(data = quants,
                aes(yintercept = cdf, colour = factor(alpha))) +
            labs(title = sprintf("%s with quantiles indicated", gp$labels$title),
               colour = "Quantile prob") +
            ## Add a cutoff to focus on smaller quantiles.
            ylim(0, quants[, max(cdf)])
    }
    gp
}

## Plot probabilities over m for each s value.
plotsprob <- function(ddist, ratio = FALSE) {
    dd <- if(ratio) {
        ddist[, list(s = s[-1], p = pmf[-1]/pmf[-.N]), by = m]
    } else {
        ddist[, list(s, p = pmf, m)]
    }
    ggplot(dd, aes(x = m, y = p)) + geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ s, scales = "free_y") +
        theme(text = element_text(size = 6)) +
        labs(title = sprintf("P[A(m) = s]%s for each s = 0:%s",
            if(ratio) "/P[A(m) = s-1]" else "", ddist[, max(m)]),
            x = "m", y = sprintf("Prob%s", if(ratio) " ratio" else ""))
}

## Plot cumulative probabilities over m for each s value.
plotscumprob <- function(ddist, ratio = FALSE) {
    dd <- if(ratio) {
        ddist[, list(s = s[-1], p = cdf[-1]/cdf[-.N]), by = m]
    } else {
        ddist[, list(s, p = cdf, m)]
    }
    ggplot(dd, aes(x = m, y = p)) + geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ s, scales = "free_y") +
        theme(text = element_text(size = 6)) +
        labs(title = sprintf("P[A(m) <= s]%s for each s = 0:%s",
            if(ratio) "/P[A(m) <= s-1]" else "", ddist[, max(m)]),
            x = "m", y = sprintf("Prob%s", if(ratio) " ratio" else ""))
}

## Plot sequence of quantiles over m for each alpha value.
plotquants <- function(quants, together = TRUE) {
    dd <- rbind(quants[, list(m, alpha, s = qs, type = "exact")],
        quants[, list(m, alpha, s = qs.interp, type = "interpolated")])
    aesvals <- list(x = quote(m), y = quote(s))
    if(together) aesvals$colour <- quote(factor(alpha))
    gp <- ggplot(dd, aes_(aesvals)) +
        geom_line(aes(linetype = type)) +
        geom_point(size = 0.5)
    if(together) gp <- gp + facet_wrap(~ alpha)
    gp + theme(text = element_text(size = 8)) +
        labs(title = "Quantile q(m; alpha) for each alpha",
            x = "m", y = "Quantile value (s)", colour = "alpha",
            linetype = "Type")
}

## Plot sequence of probabilities at each quantile over m, for each alpha.
plotquantprobs <- function(quants, exact = TRUE) {
    dd <- rbind(quants[, list(m, alpha, p = pmf, type = "exact")],
        quants[, list(m, alpha, p = pmf.interp, type = "interpolated")])
    ggplot(if(exact) dd else dd[type == "interpolated"],
           aes(x = m, y = p)) +
        geom_line(aes(linetype = type)) +
        geom_point(size = 0.5) +
        facet_wrap(~ alpha, scales = "free_y") +
        theme(text = element_text(size = 8)) +
        labs(title = "Probability at each quantile q(m; alpha) for each alpha",
            x = "m", y = "P[A(m) = q(m; alpha)]", linetype = "Type")
}

## Plot sequence of cumulative probabilities at each quantile over m,
## for each alpha.
plotquantcumprobs <- function(quants) {
    dd <- rbind(quants[, list(m, alpha, p = cdf, type = "exact")],
        quants[, list(m, alpha, p = alpha, type = "interpolated")])
    ggplot(dd, aes(x = m, y = p)) +
        geom_line(aes(linetype = type)) +
        geom_point(size = 0.5) +
        facet_wrap(~ alpha) +
        theme(text = element_text(size = 8)) +
        labs(title = paste("Cumulative probability at each quantile q(m; alpha)",
                "for each alpha"),
            x = "m", y = "P[A(m) <= q(m; alpha)]", linetype = "Type")
}

## Plot the sequence of privacy ratio values at each quantile over m,
## for each alpha.
plotquantpr <- function(prquants, exact = TRUE) {
    ## Add indicator of whether previous exact quantile was same.
    prquants[, sameasprev := c(FALSE, qs[-1] == qs[-.N]), by = alpha]
    dd <- rbind(prquants[, list(m, alpha, pr = pr, sameasprev, type = "exact")],
        prquants[, list(m, alpha, pr = pr.interp, sameasprev = FALSE,
            type = "interpolated")])
    dd[, pointsize := ifelse(type == "exact", 1, 0.5)]
    ggplot(if(exact) dd else dd[type == "interpolated"],
            aes(x = m, y = pr)) +
        geom_line(aes(linetype = type)) +
        geom_point(aes(colour = sameasprev, size = type)) +
        scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        scale_size_manual(values = c("exact" = 1, "interpolated" = 0.5),
            guide = FALSE) +
        facet_wrap(~ alpha, scales = "free_y") +
        theme(text = element_text(size = 8)) +
        labs(title = sprintf(paste("Privacy ratio value at each quantile",
            "q(m; alpha) for each alpha (n = %s, q = %s)"), n, q),
            x = "m", y = "PR(q(m; alpha); m)", linetype = "Type",
            colour = "Q(m) = Q(m-1)?")
}

pdf(sprintf("%ss-quantiles.pdf", pvar), 12, 10); gp; dev.off()

pdf("privratio-byquantiles-n100.pdf", 12, 8)
plotquantpr(prquants, FALSE)
dev.off()

