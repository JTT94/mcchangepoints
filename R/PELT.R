#' Pruned Exact Linear Time
#' @param  data = sequence of data
#' @param  measure = a function measuring error of fiting (higher value = worse fitting)
#' @param  beta = penalty
#' @param  K = pruning penalty
#' @param  markov = whether we are modelling the data in markov, if yes, we return the trans prob of segments as well
#' 
#' @return  a list containing i) vector of change points ii) segments of data iii) transition prob of segments

PELT <- function(data, measure, beta, K, markov = TRUE){
  
  # Initialization
  n = length(data) 
  f = rep(0, n + 1); f[1] = -beta;
  cp = rep(list(c()), n + 1)
  R = rep(list(c()), n+1); R[[2]] = c(1)
  
  # Main Algorithm
  for (t_star in 2:(n)){
    # step 1
    hold = c();
    for (item in length(R[[t_star]])){
      hold = c(hold, f[R[[t_star]][item]] + measure(data[(R[[t_star]][item]+1):t_star]) + beta)
    }
    f[t_star] = min(hold)
    
    # step 2
    tau = R[[t_star]][which.min(hold)]
    
    # step 3
    cp[[t_star]] = c(cp[[tau]], tau)
    
    # step 4
    for (tau in R[[t_star]]){
      if (f[tau] + measure(data[(tau + 1) : t_star]) + K <= f[t_star]){
        R[[t_star + 1]] = c(R[[t_star + 1]], tau)
      }
    }
    R[[t_star + 1]] = c(R[[t_star + 1]], t_star)
  }
  if (length(cp[[n]]) == 1){
    print("There are no change points detected")
  }
  else{
    if (markov){
      points = cp[[n]][-1]
      ls = rep(list(c()), length(cp[[n]]) - 1)
      P_ls = rep(list(), length(cp[[n]] - 1))
      cp = c(1, cp[[n]][-1], length(data))
      for (i in 1: (length(cp) - 1)){
        ls[[i]] = (data[(cp[i]+1):cp[i+1]])
        P_ls[[i]] = markovchainFit(data[(cp[i]+1):cp[i+1]])$estimate
      }
      return(list(cp = points, seg = ls, P_ls = P_ls))
    }
    else{
      points = cp[[n]][-1]
      ls = rep(list(c()), length(cp[[n]]) - 1)
      cp = c(1, cp[[n]][-1], length(data))
      for (i in 1: (length(cp) - 1)){
        ls[[i]] = (data[(cp[i]+1):cp[i+1]])
      }
      return(list(cp = points, seg = ls, P_ls = NULL))
    }
  }
}
