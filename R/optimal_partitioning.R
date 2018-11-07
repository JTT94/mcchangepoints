#' Optimal Partitioning
#' @param  data = sequence of data
#' @param  measure = a function measuring error of fiting (higher value = worse fitting)
#' @param  beta = penalty
#' @param  markov = whether we are modelling the data in markov, if yes, we return the trans prob of segments as well
#'
#' @return  a list containing i) vector of change points ii) segments of data iii) transition prob of segments

opt = function(data, measure, beta, markov = TRUE){
  
  # Initialization
  n = length(data)
  f = rep(0, n+1); f[1] = -beta # f is the optimization objective
  cp = rep(list(c()), n + 1)
  
  # Dynamic Programming
  for (t_star in 2:n){
    hold = c()
    for (i in 1:(t_star - 1)){
      hold = c(hold, f[i] + measure(data[(i + 1):t_star]) + beta)
    }
    f[t_star] = min(hold)
    tau = which.min(hold)
    cp[[t_star]] = c(cp[[tau]], tau)
  }
  
  # Results
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
