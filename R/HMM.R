
#==============================================================================
# Baum-Welch algorithm
#==============================================================================

#------------------------------------------------------------------------------

# Function to add 2 numbers on the log scale
# Uses: log(a+b) = log(a) + log(1 + b/a)
# Assumes 'a' and 'b' are log-scaled quantities
log_add <- function(a, b) {
  if(a > b) {
    ans <- a + log(1 + exp(b - a))
  }
  else {
    ans <- b + log(1 + exp(a - b))
  }
  return(ans)
}


# Implements the Forward algorithm
# All calculations are done on a log scale as probabilities decay to near 0 very quickly
log_forward_alg <- function(d, pi, A, B) {
  nStates <- dim(A)[1]
  nData <- length(d)

  # Initialising f for t=1
  f <- array(NA, c(nStates, nData))
  f[,1] <- log(pi) + log(B[d[1], d[1],])

  for (t in 2:nData) {
    for (i in 1:nStates) {
      logsum <- -Inf

      for(j in 1:nStates) {
        # Calculate next term in the summation
        logterm <- f[j,t-1] + log(A[j,i])

        # If the next term is nonzero, add it to the sum
        if(logterm > -Inf) {
          logsum <- log_add(logsum, logterm)
        }
      }
      # Calculate f[i,t] = B(d[t-1], d[t], i) * SUM_j(f[j, t-1] * A[j,i])
      f[i,t] <- log(B[d[t-1],d[t],i]) + logsum
    }
  }
  f
}

# Implements the Backward algorithm
log_backward_alg <- function(d, A, B) {
  nStates <- dim(A)[1]
  nData <- length(d)

  b <- array(NA, c(nStates, nData))
  b[,nData] <- 0

  # This loops  over the data and the states and calculates the sum in log scale
  for (t in (nData-1):1) {
    for (i in 1:nStates) {
      logsum = -Inf
      for(j in 1:nStates) {
        logterm <- b[j, t+1] + log(A[i,j]) + log(B[d[t], d[t+1],j])
        if (logterm > -Inf) {
          logsum <- log_add(logsum, logterm)
        }
      }
      b[i,t] <- logsum
    }
  }
  b
}

# Calculates P(data) as sum_i f_i(T)
log_denom <- function(f) {
  nStates <- dim(f)[1]
  nData <- dim(f)[2]

  logsum <- f[1,nData]

  for(i in 2:nStates) {
    logterm <- f[i, nData]
    if(logterm > -Inf) {
      logsum <- log_add(logsum, logterm)
    }
  }
  logsum
}

# Another way to calculate P(data) (for checking consistency) as sum_i f_i(t)*b_i(t) for any t
log_denom_2 <- function(f,b) {
  nStates <- dim(f)[1]
  nData <- dim(f)[2]

  probData <- array(NA, nData)

  for(t in 1:nData) {
    logsum <- f[1, t] + b[1,t]
    for(i in 2:nStates) {
      logterm <- f[i, t] + b[i,t]
      if(logterm > -Inf) {
        logsum <- log_add(logsum, logterm)
      }
    }
    probData[t] <- logsum
  }
  probData
}


# Calculates the posterior probability of being in state i at time t
log_gamma <- function(f, b, logdenom) {
  lg <- array(NA, dim(f))
  lg <- f + b - logdenom
}

# Calculates update of A according to the Baum-Welch alg
A_update <- function(d, f, b, A, B, logdenom) {
  nStates <- dim(f)[1]
  nData <- dim(f)[2]

  A_new <- array(0, dim(A))
  A_new[nStates, nStates] <- 1

  for (i in 1:(nStates-1)) {
    for (j in i:(i+1)) {
      logsum <- f[i,1] + log(A[i,j]) + b[j,2] + log(B[d[1],d[2],j])
      for(t in 2:(nData-1)) {
        logterm <- f[i,t] + log(A[i,j]) + b[j,t+1] + log(B[d[t],d[t+1],j])
        if(logterm > -Inf) {
          logsum <- log_add(logsum, logterm)
        }
      }
      A_new[i,j] <- exp(logsum - logdenom)
    }
  }
  A_new
}

# Calculate update of B according to the Baum-Welch alg
B_update <- function(d, f, b, A, B, logdenom, verbose = FALSE) {
  nStates <- dim(f)[1]
  nData <- dim(f)[2]
  nObs <- dim(B)[1]
  B_new <- array(NA, dim(B))

  for (j in 1:nStates) {
    for (obs in 1:nObs) {
      for (prevobs in 1:nObs) {
      logsum <- -Inf
      for (t in 2:nData) {
        if (obs == d[t] & prevobs == d[t-1]) {
          logterm <- f[j, t] + b[j, t]
          if(logterm > -Inf) {
            logsum <- log_add(logsum, logterm)
          }
        }
      }
      if (verbose) print(paste0("prevobs ", prevobs, ", obs ", obs, ", j ", j))
      if (verbose) print(exp(logsum-logdenom))
      B_new[prevobs, obs, j] <- exp(logsum - logdenom)
      }
    }
  }
  B_new
}

# One step of updates with the B-W algorithm
BW_step <- function(d, pi, A, B) {
  f <- log_forward_alg(d, pi, A, B)
  b <- log_backward_alg(d, A, B)
  logdenom <- log_denom(f)

  A_new <- A_update(d, f, b, A, B, logdenom)
  B_new <- B_update(d, f, b, A, B, logdenom)
  A_new <- (A_new/apply(A_new, 1, sum))
  for (j in 1:dim(B)[3]) {
    B_new[,,j] <- (B_new[,,j]/apply(B_new[,,j], 1, sum))
  }

  pi_new <- exp(log_gamma(f,b,logdenom)[,1])

  return(list(A_new = A_new, B_new = B_new, pi_new = pi_new))
}

# B-W algorithm recursion
BW_recursion <- function(d, pi, A, B, max_iter=100, delta = 1e-06) {

  A_temp <- A
  B_temp <- B
  pi_temp <- pi
  distances <- c()
  ll <- c()

  f <- log_forward_alg(d, pi_temp, A_temp, B_temp)
  ll[1] <- log_denom(f)

  for (i in 1:max_iter) {
    print(i)

    bw <- BW_step(d, pi_temp, A_temp, B_temp)
    A_new <- bw$A_new
    B_new <- bw$B_new
    pi_new <- bw$pi_new

    distance <- sqrt(sum((A_temp-A_new)^2)) + sqrt(sum((B_temp-B_new)^2))
    distances <- c(distances, distance)

    f <- log_forward_alg(d, pi_new, A_new, B_new)
    ll[i+1] <- log_denom(f)

    A_temp <- A_new
    B_temp <- B_new
    pi_temp <- pi_new

    if(distance < delta) {
      break
    }

  }
  return(list(A_update = A_temp, B_update = B_temp, pi_update = pi_temp, distances = distances, ll = ll))
}

#==============================================================================
# Viterbi algorithm
#==============================================================================

# This follows the same structure as the 'viterbi' function given in the HMM package:
#         https://cran.r-project.org/web/packages/HMM/HMM.pdf
# but modifies it slightly to allow B to be 'nObs x nObs x nStates'-dimensional

# Finds most likely path given data, pi, A, B
viterbi <- function (d, pi, A, B) {
  nObs <- length(d)
  nStates <- dim(A)[1]
  v <- array(NA, c(nStates, nObs))

  for (state in 1:nStates) {
    v[state, 1] = log(pi[state])
  }

  for (t in 2:nObs) {
    for (state in 1:nStates) {
      maximum <- NULL
      for(prevstate in 1:nStates) {
        temp <- v[prevstate, t-1] + log(A[prevstate, state])
        maximum <- max(maximum, temp)
      }
      v[state, t] <- log(B[d[t-1],d[t],state]) + maximum
    }
  }

  path <- rep(NA, nObs)
  for(state in 1:nStates) {
    if(max(v[, nObs]) == v[state, nObs]) {
      path[nObs] <- state
      break
    }
  }

  for (t in (nObs-1):1) {
    for (state in 1:nStates) {
      if(max(v[, t] + log(A[, path[t+1]])) == v[state, t] + log(A[state, path[t+1]])) {
        path[t] <- state
        break
      }
    }
  }

  path
}







