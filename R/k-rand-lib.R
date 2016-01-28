#############################################################################
###
###  Useful functions for running simulations for k-randomization project.
###
#############################################################################


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

## Compute the privacy ratio for the case L = 1:
## - for a collection of size n and lie probability q
## - where the original collection has m 1s
## - and the modification changes a 1 to a 0.
## The privacy ratio P[A(m-1) = s]/P[A(m) = s] is computed across all
## outcome values s = 0,...,n and over all original collections m = 1,...,n.
## As in the notation of section 3.1, A is a rv with distribution 
## Bin(m, 1-q) + Bin(N-m, q).
## Returns a list of length n indexed by m, whose elements are vectors of
## length n+1 giving the privacy ratio values for that m indexed by s+1.
privratio <- function(n, q) {
    ## Compute the full pmf of A(m) for all values of m.
    pmfA <- lapply(0:n, dbinomsum, n, q)
    ## The privacy ratio is pmfA[[m]]/pmfA[[m+1]] for m = 1,...,n.
    lapply(1:n, function(m) { pmfA[[m]] / pmfA[[m+1]] })
}


## Find the alpha-th quantile of the distribution of Y as in dbinomsum(m,n,q).
## This is the largest number x such that P[Y < x] <= alpha, ie.
## P[Y >= x] >= 1 - alpha.
qbinomsum <- function(m, n, q, alpha) {
    ## First compute the cdf.
    cdf <- cumsum(dbinomsum(m, n, q))
    ## Find the largest x such that alpha >= P[Y <= x-1] = P[Y < x].
    lta <- cdf <= alpha
    if(!any(lta)) return(0)
    max(which(cdf <= alpha))
}

