
fit_single_matrix <- function(chain) {
  attr(markovchain::markovchainFit(chain)$estimate, "transitionMatrix")
}

fit_multiple_matrices <- function(chain, change_points) {

  # Check change points
  if (!length(change_points)>0) {
    return( list(fit_single_matrix(chain) ))
  }
  change_points <- sort(change_points)
  if (! head(change_points, 1) %in% c(0)) {
    change_points <- c(0,change_points)
  }

  if (! tail(change_points, 1) >= length(chain)) {
    change_points <- c(change_points, length(chain))
  }

  # get chain index
  index <- seq_along(chain)

  #cut up index by change points buckets
  cut_index <- cut(index, change_points, labels = F)

  # fit matrices by cut index
  mle_mats <- list()
  for (i in unique(cut_index)) {
    chain_subset <- chain[cut_index == i]
    mle_mats[[i]] <- fit_single_matrix(chain_subset)
  }

  mle_mats
}


