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

L <- 2
ntypes <- 2^L
## Try computing distributions for L=2.
## Assume the multinomial categories 1, 2, 3, 4 correspond to
## (11), (10), (01), (00) respectively.

## Take a given initial configuration.
initconfig <- c(3, 2, 3, 2)
n <- sum(initconfig)
## Create a table of all possible synthetic collections.
synthcoll <- as.data.table(do.call(expand.grid, rep(list(1:ntypes), n)))
synthcollcols <- sprintf("x%s", 1:n)
setnames(synthcoll, synthcollcols)
## Summarize collections by mapping to multinomial counts.
bincols <- sprintf("inbin%s", 1:n)
mncols <- sprintf("s%s", 1:ntypes)
for(j in 1:ntypes) {
    synthcoll[, eval(bincols) := lapply(synthcollcols, function(v) { 
        get(v) == j })]
    synthcoll[, eval(mncols[[j]]) := Reduce("+", lapply(bincols, function(bc) {
        synthcoll[[bc]] }))]
    synthcoll[, eval(bincols) := NULL]
}
## Now compute probabilities of obtaining each synthetic collection.
## First find the distance metric between the original and synthetic.
distvals <- list(c(0, 1, 1, 2), c(1, 0, 2, 1), c(1, 2, 0, 1), c(2, 1, 1, 0))
distfn <- function(v, j) { distvals[[j]][v] }
initcoll <- rep(1:4, each = initconfig)
synthcoll[, d := {
    dvals <- mapply(function(origv, scol) { distfn(get(scol), origv) },
        initcoll, synthcollcols, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    Reduce("+", dvals) }]
synthcoll[, p := (1-q)^(n*L)*(q/(1-q))^d, by = d]

## The probability of getting any multinomial outcome can now be obtained by
## summing the relevant rows in the table.
mnprobs <- synthcoll[, list(p = sum(p)), by = mncols]
setkeyv(mnprobs, mncols)

## Compute the probability ratio for a given sequence of s values, given as
## a data table with 4 columns named by mncols.
probratio <- function(s, i, j) {
    num <- mnprobs[s, p]
    sdenom <- copy(s)[, eval(mncols[[i]]) := get(mncols[[i]]) - 1][,
       eval(mncols[[j]]) := get(mncols[[j]]) + 1]
    denom <- mnprobs[sdenom, p]
    num / denom
}

## Get some data.
#s <- mnprobs[s4 == 0][, p := NULL]
s <- mnprobs[s3 == 0][, p := NULL]
s[, pr := probratio(s, 4, 2)]



#-------------

## Create a table of all possible multinomial outcomes.
n <- 10
mnvals <- as.data.table(do.call(expand.grid, rep(list(0:n), ntypes- 1)))
cols <- sprintf("s%s", 1:(ntypes))
setnames(mnvals, cols[-length(cols)])
mnvals[, eval(cols[ntypes]) := n - rowSums(mnvals)]
mnvals <- mnvals[get(cols[ntypes]) >= 0]

## Probability of getting a given outcome (s_1,...,s_{2^L}):
## Product of 2^L multinomial probabilities, summed over all subsets of size
## 2^L of rows of mnvals such that the j-th column sum is equal to s_j.

## Figure out general form for probabilities later on.
## Matrix of probabilities for each component multinomial.
mnprobs <- rbind(
    c((1-q)^2, q*(1-q), q*(1-q), q^2),
    c(q*(1-q), (1-q)^2, q^2, q*(1-q)),
    c(q*(1-q), q^2, (1-q)^2, q*(1-q)),
    c(q^2, q*(1-q), q*(1-q), (1-q)^2))

## Compute the probability of generating a vector with values (s_1,...,s_{2^L})
## from a sum of 2^L multinomials, with m_i trials each.
dmnsum <- function(s, m, mnp) {
    ## Find all combinations outcomes for each multinomial that together
    ## will combine to give the required s.
    suboutcomes <- lapply(s, function(v) {
        vv <- as.data.table(do.call(expand.grid, rep(list(0:v), ntypes)))
        vv[rowSums(vv) == v]
    })
}


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

