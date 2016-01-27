#############################################################################
###
###  Useful functions for running simulations for k-randomization.
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
    ## In all probabilities below, P[Bin(K,r) = j] = 0 if j < 0 or j > K.
    ## [P[Bin(m,p) = j] for j = 0,...,n].
    bin.mp <- dbinom(0:n, m, 1-q)
    ## [P[Bin(n-m,q) = j] for j = 0,...,n].
    bin.nmq <- dbinom(0:n, n-m, q)
    ## Create a (n+1) x (n+1) matrix where row i is
    ## [P[Bin(n-m,q) = i], P[Bin(n-m,q) = i-1], ..., P[Bin(n-m,q) = i-n]].
    convol <- shift(rev(bin.nmq), n = n:0, fill = 0, type = "lead")
    convol <- do.call(rbind, convol)
    ## The result of the convolution is a vector where entry i is
    ## sum_{j=0}^n P[Bin(n-m,q) = i-j]*P[Bin(m,p) = j].
    as.vector(convol %*% bin.mp)
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

