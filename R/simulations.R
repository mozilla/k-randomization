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

## Computations for the general case.

## Bitwise randomization probability.
q <- 0.2
## The dimesion of the bit vectors.
L <- 2
## The number of unique vectors in the space.
ntypes <- 2^L
## Vector-wise transition probabilities for applying the randomization.
pmat <- transmatrix(L, q)

## The configuration of the original collection.
m <- c(10, 5, 6, 9)
## Compute the distribution P[A(m) = s].
mndist <- mnprobs(m, pmat)
## Column names, for reference.
scols <- scolnames(mndist)
newscols <- sprintf("new%s", scols)

i <- 1
j <- 4
pr <- probratio(mndist, i, j)


## For a given e_{kr} shift vector, compute side-by-side probability ratio
## values for the original and shifted vectors.
## Returns a table containing columns "prorig" and "prnew".
## Does not modify the original pr table.
compareprs <- function(pr, from, to) {
    prcomp <- copy(pr)
    setnames(prcomp, "pr", "prorig")
    ## Compute the shifted s values.
    shiftsvals(prcomp, from, to)
    ## Keep only valid ones.
    setkeyv(prcomp, newscols)
    prcomp <- prcomp[pr, nomatch = 0]
    setnames(prcomp, "pr", "prnew")
    prcomp
}

## Check monotonicity and maximality of the probability ratio from i to j in all
## relevant directions, for all valid s values:
## - Moving away from i is increasing
##      > rho(s, i, j) <= rho(s + e_{ik}, i, j)
## - Moving towards j when i is constant is increasing
##      > rho(s, i, j) <= rho(s + e_{kj}, i, j), k != j
## - Moving from i to j increases more than moving away from i in another
##   direction
##      > rho(s + e_{ik}, i, j) <= rho(s + e_{ij}, i, j), k != j
## - Moving from i to j increases more than moving towards j from another
##   direction when i is constant
##      > rho(s + e_{kj}, i, j) <= rho(s + e_{ij}, i, j), k != j
##
checkprproperties <- function(pr, i, j) {
    cat(sprintf("Parameters: 2^L = %s, n = %s, m = %s, i = %s, j = %s\n",
        ntypes, n, sprintf("(%s)", paste(m, collapse = ",")), i, j))
    k <- 1:ntypes
    k <- k[k != i]
    awayfromi <- as.logical(lapply(k, function(ind) {
        compareprs(pr, i, ind)[, all(prnew >= prorig)]
    }))
    if(all(awayfromi)) cat("Moving away from i always increasing\n")
    else cat(sprintf("Moving away from i not increasing when k = %s\n",
        paste(k[!awayfromi], collapse = ",")))

    ## k is now all indices except i or j.
    k <- k[k != j]
    towardsj <- as.logical(lapply(k, function(ind) {
        compareprs(pr, ind, j)[, all(prnew >= prorig)]
    }))
    if(all(towardsj)) cat("Moving towards j always increasing\n")
    else cat(sprintf("Moving towards j not increasing when k = %s\n",
        paste(k[!towardsi], collapse = ",")))

    ## Precompute the probability ratio on shifting from i to j.
    prij <- compareprs(pr, i, j)
    setnames(prij, "prnew", "prij")
    setkeyv(prij, scolnames(prij))
    maxawayfromi <- as.logical(lapply(k, function(ind) {
        prcomp <- compareprs(pr, i, ind)
        setkeyv(prcomp, scolnames(prcomp))
        prij[prcomp, nomatch = 0][, all(prij >= prnew)]
    }))
    if(all(maxawayfromi))
        cat("Moving away from i increases most in j direction\n")
    else cat(sprintf("Moving away from i increases more when k = %s than j\n",
        paste(k[!maxawayfromi], collapse = ",")))

    maxtowardsj <- as.logical(lapply(k, function(ind) {
        prcomp <- compareprs(pr, ind, j)
        setkeyv(prcomp, scolnames(prcomp))
        prij[prcomp, nomatch = 0][, all(prij >= prnew)]
    }))
    if(all(maxtowardsj))
        cat("Moving towards j increases most from j direction\n")
    else cat(sprintf("Moving towards j increases more when k = %s than i\n",
        paste(k[!maxtowardsj], collapse = ",")))
}

#----------------------------------------

shiftv <- function(k,l) {
    v <- rep(0, ntypes)
    v[k] <- -1
    v[l] <- 1
    v
}

## Compare the probability ratio computed from mndist to that computed using the
## recursion, conditioning on an original vector of type r.
prr <- function(s, i, j) { prrecursion(s, i, j, 1, m, pmat, verbose = TRUE) }




##########################################################

## Graphs for L=1 case.

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


## Compute pmfs for m = 0,...,n.
pmfs <- allpmfs(n, q)
## Tabulate distribution (pmf/cdf) over m,s.
## Creates a table of (m, s, pmf, cdf).
ddist <- rbindlist(mapply(function(i, v) {
    data.table(m = i - 1, s = seq_along(v) - 1, pmf = v)
}, seq_along(pmfs), pmfs, SIMPLIFY = FALSE, USE.NAMES = FALSE))
ddist[, cdf := cumsum(pmf), by = m]
## Add shifted versions of CDFs for each m.
ddist[, ssf := s - floor(m*(1-2*q)), by = m]
ddist[, ssc := s - ceiling(m*(1-2*q)), by = m]



## Compute privacy ratio for all s, m = 0,...,n.
pr <- privratio(pmfs)
## Compute table of quantile info (distribution values at quantiles).
quants <- pmfquantile(pmfs, alpha, interpolate = TRUE)

## Privacy ratio at quantiles.
setkey(pr, m, s)
setkey(quants, m, qs)
prquants <- quants[pr][, list(m, alpha, qs, pr, qs.interp, interp)]
prquants[, pr.interval := rev(diff(c(rev(pr), 0))), by = m]
prquants[, pr.interp := pr + interp * pr.interval][, pr.interval := NULL]
prquants <- prquants[!is.na(alpha)]


## Difference between shifted distributions, compared to CDF jump sizes.
ddiff <- ddist[s <= max(ssf[m == n]) & m == 0][, list(s, cdf0 = cdf)]
setkey(ddiff, s)
## Compare to F_0 shifted up or down by 1 or 2.
ddiff[, lcdf0 := c(0, cdf0[-length(cdf0)])][,
    l2cdf0 := c(0, 0, cdf0[-(length(cdf0) + (-1:0))])][,
    ucdf0 := c(cdf0[-1], 1)][,
    u2cdf0 := c(cdf0[-(1:2)], 1, 1)]
## Consecutive differences for F_0.
#ddiff[, ldiff := c(cdf[1], diff(cdf))][,
#    ldiff2 := c(cdf[1:2], diff(cdf, lag = 2))][,
#    udiff := c(diff(cdf), 1-cdf[length(cdf)])][,
#    udiff2 := c(diff(cdf, lag = 2), 1-c(tail(cdf, 2)))]
#setnames(ddiff, "cdf", "cdf0")
## Compare against shifted versions of F_m.
srange <- ddiff[, range(s)]
fshift <- ddist[m > 0 & ssf >= srange[1] & ssf <= srange[2]][,
    list(m, s = ssf, fcdf = cdf)]
setkey(fshift, s)
## Join against table for F_0.
ddiff <- ddiff[fshift]
## Same for upper shift, but join on both s and m.
setkey(ddiff, s, m)
cshift <- ddist[m > 0 & ssc >= srange[1] & ssc <= srange[2]][,
    list(m, s = ssc, ccdf = cdf)]
setkey(cshift, s, m)
ddiff <- ddiff[cshift]

## Try same thing for PMFs.
ddiff <- ddist[s <= max(ssf[m == n]) & m == 0][, list(s, pmf0 = pmf)]
setkey(ddiff, s)
ddiff[, lpmf0 := c(0, pmf0[-length(pmf0)])]
## Compare against shifted versions of p_m.
srange <- ddiff[, range(s)]
fshift <- ddist[m > 0 & ssf >= srange[1] & ssf <= srange[2]][,
    list(m, s = ssf, fpmf = pmf)]
setkey(fshift, s)
## Join against table for F_0.
ddiff <- ddiff[fshift]
## Same for upper shift, but join on both s and m.
setkey(ddiff, s, m)
cshift <- ddist[m > 0 & ssc >= srange[1] & ssc <= srange[2]][,
    list(m, s = ssc, ccdf = cdf)]
setkey(cshift, s, m)
ddiff <- ddiff[cshift]

## Look at ratios of consecutive PMF values.
probr <- ddist[, list(s = s[-1], pr = pmf[-1]/pmf[-(n+1)]), by = m]
mum <- data.table(m = 0:n)[, mu := n*q + m*(1-2*q)]
deltas <- mum[m == 0, 0:floor(mu - 1)]
mus <- mum[, list(delta = deltas, s = floor(mu - deltas)), by = m]
setkey(mus, m, s)
setkey(probr, m, s)
probrm <- probr[mus]


##########################################################

## Plot CDFs for multiple m.

## Plot cdfs for multiple m on the same graph.
## If shift = TRUE, plots F(x + floor(m(p-q))).
## If additionally use.ceiling = TRUE, replaces 'floor' with 'ceiling'.
## Optionally limit the y axis range using 'lower' and 'upper'.
plotcdfshift <- function(ddist, mvals, shift = FALSE, use.ceiling = FALSE,
                                                lower = NULL, upper = NULL) {
    xx <- if(shift) {
        ## Find the maximum s range to consider over our choice of m.
        svals <- ddist[m %in% mvals,
            do.call(seq, as.list(range(if(use.ceiling) ssc else ssf)))]
        smvals <- expand.grid(svals, mvals)
        smvals <- data.table(m = smvals[[2]], s = smvals[[1]], key = c("m","s"))
        setkeyv(ddist, c("m", if(use.ceiling) "ssc" else "ssf"))
        ddshifted <- ddist[smvals]
        setnames(ddshifted, if(use.ceiling) "ssc" else "ssf", "ss")
        ddshifted[, s := NULL][is.na(pmf), pmf := 0][
            is.na(cdf) & ss < 0, cdf := 0][
            is.na(cdf), cdf := 1]
        setnames(ddshifted, "ss", "s")
        ddshifted
    } else ddist[m %in% mvals]
    ylims <- c(0,1)
    if(!is.null(lower)) ylims[1] <- lower
    if(!is.null(upper)) ylims[2] <- upper
    qp <- qplot(s, cdf, colour = factor(m), data = xx, geom = "step") +
        geom_hline(aes(yintercept = cdf), data = xx[m == min(m)], 
            colour = "grey30", linetype = "dotted") +
        scale_y_continuous(breaks = interval.breaks(0.1), limits = ylims)
    titlestr <- sprintf("CDF of A(m) %%s[n = %s, q = %s]", n, q)
    if(shift) {
        qp <- qp + xlim(-5, 25)
        titlestr <- sprintf(titlestr, sprintf("shifted by %s(m(p-q))",
            if(use.ceiling) "ceiling" else "floor"))
    } else {
        titlestr <- sprintf(titlestr, "")
    }
    qp + labs(title = titlestr, colour = "m") +
        theme(text = element_text(size = 8))
}


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

## Plot the probability ratio P[s]/P[s-1] over s = 1,...,n, for each m.
plotprobratio <- function(ddist) {
    dd <- ddist[, list(s = s[-1], p = pmf[-1]/pmf[-.N]), by = m]
    ggplot(dd, aes(x = s, y = p)) + geom_line() +
        geom_point(size = 0.5) +
        facet_wrap(~ m
            #, scales = "free_y"
            ) +
        theme(text = element_text(size = 6)) +
        labs(title = sprintf("P[A(m) = s]/P[A(m) = s-1] for each s = 0:%s",
            ddist[, max(m)]),
            x = "s", y = "Prob ratio")
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
    aesvals <- list(x = "m", y = "s")
    if(together) aesvals$colour <- "factor(alpha)"
    gp <- ggplot(dd, do.call("aes_string", aesvals)) +
        geom_line(aes(linetype = type)) +
        geom_point(size = 0.5)
    if(!together) gp <- gp + facet_wrap(~ alpha)
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
plotquantpr <- function(prquants, exact = TRUE, fixedy = FALSE) {
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
        (if(fixedy) facet_wrap(~ alpha)
            else facet_wrap(~ alpha, scales = "free_y")) +
        theme(text = element_text(size = 8)) +
        labs(title = sprintf(paste("Privacy ratio value at each quantile",
            "q(m; alpha) for each alpha (n = %s, q = %s)"), n, q),
            x = "m", y = "PR(q(m; alpha); m)", linetype = "Type",
            colour = "Q(m) = Q(m-1)?")
}

## Plot ratios of ratios.
plotquantprratios <- function(prquants, absval = FALSE) {
    ## Add indicator of whether previous exact quantile was same.
    prquants[, sameasprev := c(FALSE, qs[-1] == qs[-.N]), by = alpha]
    dd <- prquants[, list(m = m[-1], rpr = pr[-1]/pr[-.N],
        sameasprev = sameasprev[-1]), by = alpha]
    if(absval) dd[sameasprev == FALSE, rpr := 1/rpr]
    ## Find max for both cases.
    ddmax <- dd[, .SD[which.max(rpr)], by = list(alpha, sameasprev)]
    gp <- ggplot(dd, aes(x = m, y = log(rpr))) + geom_line() +
        geom_point(aes(colour = sameasprev)) +
        facet_wrap(~ alpha) +
        geom_point(data = ddmax, size = 0.5) +
        theme(text = element_text(size = 8)) +
        labs(title = sprintf(paste("Log ratio of consecutive privacy ratio",
            "values at quantiles q(m; alpha) for each alpha (n = %s, q = %s)"),
                n, q),
            x = "m", y = "log(PR(q(m+1; alpha); m+1)/PR(q(m+1; alpha); m+1))",
            colour = "Q(m) = Q(m-1)?")
    if(!absval)  gp <- gp + geom_hline(yintercept = 0, colour = "grey30")
    gp
}

pdf(sprintf("%ss-quantiles.pdf", pvar), 12, 10); gp; dev.off()

pdf("privratio-byquantiles-n100.pdf", 12, 8)
plotquantpr(prquants, FALSE)
dev.off()


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
    if(!is.null(quants)) {
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

