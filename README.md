MultiBDPureBirth
======

`MultiBDPureBirth` is a slight extention of the `MultiBD`  `R` package for direct likelihood-based inference of multivariate birth-death processes. 

## Installation

1. Install (if necessary) package dependencies and helpers:
```{r}
install.packages(c("Rcpp", "RcppParallel", "BH", "devtools"))
```

2. Install `MultiBDPureBirth` from `github` (until package becomes available via `CRAN`):
```{r}
devtools::install_github("msuchard/MultiBDPureBirth")
```

## Short example

```{r}
library(MultiBDPureBirth)
data(Eyam)

loglik_sir <- function(param, data) {
  alpha <- exp(param[1]) # Rates must be non-negative
  beta  <- exp(param[2])
  
  # Set-up SIR model
  drates1 <- function(a, b) { 0 }
  brates2 <- function(a, b) { 0 }
  drates2 <- function(a, b) { alpha * b     }
  trans12 <- function(a, b) { beta  * a * b }
  
  sum(sapply(1:(nrow(data) - 1), # Sum across all time steps k
             function(k) {
               log(
                 dbd_prob(  # Compute the transition probability matrix
                   t  = data$time[k + 1] - data$time[k], # Time increment
                   a0 = data$S[k], b0 = data$I[k],       # From: S(t_k), I(t_k)                                      
                   drates1, brates2, drates2, trans12,
                   a = data$S[k + 1], B = data$S[k] + data$I[k] - data$S[k + 1],
                   computeMode = 4, nblocks = 80         # Compute using 4 threads
                 )[1, data$I[k + 1] + 1]                 # To: S(t_(k+1)), I(t_(k+1))
               )
             }))
}

loglik_sir(log(c(3.204, 0.019)), Eyam) # Evaluate at mode
```


## License
`MultiBDPureBirth` is licensed under Apache License 2.0


## References

1. Ho LST, Xu J, Crawford FW, Minin VN, Suchard MA.
[Birth(death)/birth-death processes and their computable transition probabilities with statistical applications](https://arxiv.org/abs/1603.03819).
*arXiv preprint arXiv*:1603.03819, 2016.
