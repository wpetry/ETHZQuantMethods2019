#################################################-
## Stochastic population models ----
## W.K. Petry
#################################################-
## Preliminaries ----
#################################################-
library(popdemo)

acor <- function(P){
  if(!all.equal(dim(P), c(2, 2))) stop("Autocorrelation function invalid for matrices larger than 2x2")
  if(!all(P[,1] == rev(P[,2]))) stop("Autocorrelation function invalid for asymmetric matrices")
  2 * P[1, 1] - 1
}

sim_env <- function(P, iter = 50){
  # number of possible states
  num.states <- nrow(P)
  # stores the states X_t through time
  states <- numeric(iter)
  # initialize variable for first state 
  states[1] <- 1
  for(t in 2:iter) {
    # probability vector to simulate next state X_{t+1}
    p  <- P[states[t-1], ]
    ## draw from multinomial and determine state
    states[t] <-  which(rmultinom(1, 1, p) == 1)
  }
  out <- rownames(P)[states]
  # out <- states
  return(out)
}

#################################################-
## Ex: Heads & Tails ----
#################################################-
A_t <- matrix(c(0.1, 3.0,
                0.2, 0),
              nrow = 2, byrow = TRUE)
A_h <- matrix(c(0.2, 0.2,
                1.0, 0),
              nrow = 2, byrow = TRUE)
n_0 <- c(100, 100)

# Generate a sequence of environments using a coin (tails = t, heads = h)

# [HERE]




# What is the long-term stochastic population growth rate?
n_0 <- c(1e6, 1e6)
A_list <- list(t = A_t, h = A_h)
E_iid <- matrix(0.5, nrow = 2, ncol = 2,
                dimnames = list(c("t", "h"),
                                c("t", "h")))

lams_iid <- stoch(A_list, "lambda", Aseq = E_iid, vector = n_0,
                  iterations = 2000, discard = 500)
lams_iid

eigs(A_t, "lambda")  # lambda, stable t environment
eigs(A_h, "lambda")  # lambda, stable h environment

#################################################-
## Effects of autocorrelation ----
#################################################-
## (1) Independent, identically-distributed (iid) environment
E_iid
(Esim_iid <- sim_env(E_iid, iter = 1e5))[1:50]

acor(E_iid)  # calculate autocorrelation
acf(as.numeric(as.factor(Esim_iid)), plot = FALSE)$acf[2]  # check

## (2) Positive autocorrelaton (="sticky" environments)
E_pos <- matrix(c(0.9, 0.1,
                  0.1, 0.9),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("t", "h"),
                                c("t", "h")))
E_pos
acor(E_pos)  # calculate autocorrelation
(Esim_pos <- sim_env(E_pos, iter = 1e5))[1:50]  # check
acf(as.numeric(as.factor(Esim_pos)), plot = FALSE)$acf[2]

stoch(A_list, "lambda", Aseq = E_pos, vector = n_0,
      iterations = 2000, discard = 500)

## Negative autocorrelaton (="repellent" environments)
E_neg <- matrix(c(0.1, 0.9,
                  0.9, 0.1),
                nrow = 2, byrow = TRUE,
                dimnames = list(c("t", "h"),
                                c("t", "h")))
E_neg
acor(E_neg)  # calculate autocorrelation
(Esim_neg <- sim_env(E_neg, iter = 1e5))[1:50]  # check
acf(as.numeric(as.factor(Esim_neg)), plot = FALSE)$acf[2]

stoch(A_list, "lambda", Aseq = E_neg, vector = n_0,
      iterations = 2000, discard = 500)
