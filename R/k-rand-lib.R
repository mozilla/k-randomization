#############################################################################
###
###  Useful functions for running simulations for k-randomization project.
###
#############################################################################


## Generate the distribution of Y ~ MN(m_1, p) + ... + MN(m_{2^L}, p),
## a sum of multinomials, where p is a vector of category probabilities.
## Returns a data table containing all valid multinomial outcomes over n
## trials together with their associated probabilities.
## Input is m, a vector of length 2^L whose sum is n, giving the number of
## trials from each original vector type in the original collection, and a
## transition matrix P whose i-th row gives the outcome probabilities for
## the m_i trials based on original type i vectors.
mnprobs <- function(m, P) {
    n <- sum(m)
    ntypes <- length(m)
    ## First generate a reference table of all combinations of numbers
    ## 0:n.
    ## This can be subsetted to find collections of multinomial outcomes.
    outcomes <- as.data.table(expand.grid(rep(list(0:n), ntypes)))
    setnames(outcomes, sprintf("s%s", seq_along(outcomes)))
    ## Add in row sums to be used for subsetting.
    outcomes[, rowsum := rowSums(outcomes)]
    ## The rows of interest must sum to at most 10.
    outcomes <- outcomes[rowsum <= 10]
    ## For each multinomial component, find its full pmf.
    ## Take outer product of all possible combinations of all component
    ## probabilities, then aggregate by total outcome to get the overall
    ## probabilities.
    mnpmf <- outcomes[rowsum == m[[1]]][, rowsum := NULL][, oldi := .I]
    mnpmf[, p := dmultinom(as.integer(.SD[1]), prob = P[1,]), by = oldi]
    scols <- grep("^s", names(mnpmf), value = TRUE)
    for(j in 2:ntypes) {
        ## Ignore components with no trials.
        if(m[[j]] == 0) next
        ## Compute the distribution for the next component.
        newmnpmf <- outcomes[rowsum == m[[j]]][, rowsum := NULL][, i := .I]
        newmnpmf[, p := dmultinom(as.integer(.SD[1]), prob = P[j,]), by = i]
        setnames(newmnpmf, sprintf("new%s", names(newmnpmf)))
        setkey(newmnpmf, newi)
        setkey(mnpmf, oldi)
        ## Take the cartesian product.
        combined <- CJ(oldi = mnpmf$oldi, newi = newmnpmf$newi)
        setkey(combined, oldi)
        combined <- mnpmf[combined]
        setkey(combined, newi)
        combined <- newmnpmf[combined]
        ## Find product of all combinations of probabilities.
        combined[, p := p * newp]
        ## Find the total outcomes for each category.
        combined[, eval(scols) := lapply(scols, function(ncol) {
            get(ncol) + get(sprintf("new%s", ncol)) })]
        ## Aggregate over unique outcomes.
        mnpmf <- combined[, list(p = sum(p)), by = scols][, oldi := .I]
    }
    ## The result of this should include all multinomial combinations for
    ## n trials.
    if(nrow(mnpmf) != outcomes[rowsum == 10, .N])
        stop(paste("There was an error computing the distribution for",
            "multinomial sums."))
    mnpmf[, oldi := NULL]
    mnpmf
}

## Compute the transition matrix mapping an original bit vector of length L
## to its randomized version.
## Assumes the bit vectors are ordered in decreasing order of their binary
## numeric values, eg. (11), (10), (01), (00).
## Returns a stochastic matrix with 2^L rows and columns.
transmatrix <- function(L, q) {
    probs <- c(1-q, q)
    init <- rbind(probs, rev(probs))
    if(L == 1) return(init)
    pmat <- init
    ## For L > 1, the transition matrix can be computed by repeated applying
    ## a Kronecker product with the original matrix.
    ## This corresponds to conditioning on the first bit.
    for(j in 2:L) pmat <- kronecker(pmat, init)
    return(pmat)
}

##-----------------------------------------------------------------

## Compute the distribution of Y ~ Bin(m, p) + Bin(n-m, q), where p = 1-q.
## This is interpreted as Y ~ Bin(n, q) if m = 0 and Y ~ Bin(n, p) if m = n.
## Returns a vector of length n+1 giving the pmf of Y (entry i is P[Y = i-1]).
dbinomsum <- function(m, n, q) {
    if(m == 0) return(dbinom(0:n, n, q))
    if(m == n) return(dbinom(0:n, n, 1-q))
    ## Otherwise, find the pmfs for the component densities and compute
    ## the convolution.
    ## [P[Bin(m,p) = j] for j = 0,...,n].
    bin.mp <- dbinom(0:n, m, 1-q)
    ## [P[Bin(n-m,q) = j] for j = 0,...,n].
    bin.nmq <- dbinom(0:n, n-m, q)
    ## P[Y = s] = sum_{k=max(s-(n-m),0)}^{min(s,m)} P[B(m,p)=k]P[B(n-m,q)=s-k].
    ## First compute the bounds on the sum for each s = 0:n.
    upper <- pmin(0:n, m)
    lower <- pmax(0:n-n+m, 0)
    ## Compute the convolution sum for each s.
    mapply(function(s, lb, ub) {
        k <- lb:ub
        sum(bin.mp[k+1] * bin.nmq[s-k+1])
    }, 0:n, lower, upper)
}


## Compute the pmfs of Y ~ Bin(m, p) + Bin(n-m, q) for all m in 0,...,n.
## Returns a list of length n+1, of which entry i is a vector of length
## n+1 giving the pmf for m = i-1.
allpmfs <- function(n, q) {
    lapply(0:n, dbinomsum, n, q)
}


## Compute the privacy ratio for the case L = 1:
## - for a collection of size n and lie probability q
## - where the original collection has m 1s
## - and the modification changes a 1 to a 0.
## The privacy ratio P[A(m) = s]/P[A(m+1) = s] is computed across all
## outcome values s = 0,...,n and over all original collections m = 0,...,n-1.
## As in the notation of section 3.1, A is a rv with distribution 
## Bin(m, 1-q) + Bin(N-m, q).
## Takes the output of allpmfs (a list of length n+1) for the required n
## and q, and returns a data table of (m, s, pr).
## Some of the privacy ratio values may evaluate as NaN when they are the
## ratio of two very small numbers. These are replaced with NA.
privratio <- function(pmfs) {
    ## The privacy ratio is pmfA[[m]]/pmfA[[m+1]] for m = 1,...,n.
    n <- length(pmfs) - 1
    d <- rbindlist(lapply(1:n, function(m) {
        data.table(m = m, s = 0:n, pr = pmfs[[m]] / pmfs[[m+1]])
    }))
    ## Remove NaNs.
    d[is.nan(pr), pr := NA]
    d
}


## Find the alpha-th quantiles of given pmfs.
## Each pmf should be a vector whose i-th element is P[Y = i-1].
## Accepts either a single pmf or a list as returned by allpmfs().
## For each pmf, compute the alpha-th quantile for each value of alpha.
##
## Returns a data table of (m, alpha, qs, pmf, cdf), where m indexes the pmf
## and qs is the smallest integer such that P[Y_m <= qs] >= alpha.
## The columns 'pmf' and 'cdf' give probabilities at the quantile:
## P[Y_m = qs] and P[Y_m <= qs] respectively.
##
## If 'interpolate' is TRUE, the quantile is interpolated linearly between
## the nearest integer values. Additional returned columns are 'interp' (the
## fraction (P[Y_m <= qs] - alpha)/P[Y_m = qs]), 'qs.interp' (qs - interp),
## and the pmf value at qs.interp approximated by interpolating linearly down.
## Note that this is not a real probability since the pmf is discrete.
pmfquantile <- function(pmfs, alpha, interpolate = FALSE) {
    if(!is.list(pmfs)) pmfs <- list(pmfs)
    rbindlist(mapply(function(m, v) {
        ## Compute cdf from pmf.
        cdf <- cumsum(v)
        quants <- lapply(alpha, function(a) {
            ## Find the smallest cdf element that is no smaller than a.
            qs <- min(which(cdf >= a))
            ## The quantile value is 1 less (between 0 and n).
            res <- list(m = m - 1, alpha = a, qs = qs - 1, pmf = v[qs],
                cdf = cdf[qs])
            if(interpolate) {
                lb <- if(qs == 1) 0 else cdf[qs - 1]
                interp <- (cdf[qs] - a) / (cdf[qs] - lb)
                res[["qs.interp"]] <- res[["qs"]] - interp
                res[["pmf.interp"]] <-
                    v[qs] - interp * (v[qs] - if(qs == 1) 0 else v[qs-1])
                res[["interp"]] <- interp
            }
            res
        })
        rbindlist(quants)
    }, seq_along(pmfs), pmfs, USE.NAMES = FALSE, SIMPLIFY = FALSE))
}

