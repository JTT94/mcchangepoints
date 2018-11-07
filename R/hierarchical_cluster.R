# Hierarchical Clustering


hier_markov_chain_cluster <- function(chain, threshold, initial_width = 100, dist_func = frob_norm) {
  #initial change points
  number_of_observations <- length(chain)
  cut_points <- seq(0, number_of_observations, initial_width)
  merge_bool <- T

  states <- as.character(unique(chain))
  while (any(merge_bool) & (length(cut_points) > 1)) {
    print(paste0("Length cut_points: ", length(cut_points)))
    print(paste0("Length merge_bool: ", sum(merge_bool)))
    # fit matrices
    mats <- fit_multiple_matrices(chain, change_points = cut_points)
    mats <- lapply(mats, function(x) standardize_matrix(cols = states, rows = states, mat=x))

    # merge blocks
    distance <- c()
    merge_bool <- c()
    for(i in seq(1, length(mats)-1, 2)) {
      dist <- dist_func(mats[[i]], mats[[i+1]])
      distance <- c(distance,dist)
      merge_bool <- c(merge_bool, dist < threshold)
    }
    cut_points <- cut_points[!merge_bool]
  }
  mats <- fit_multiple_matrices(chain, change_points = cut_points)
  mats <- lapply(mats, function(x) standardize_matrix(cols = states, rows = states, mat=x))

  return( list('matrices' = mats,
               'change_points' = cut_points))
}




one_step_markov_chain_cluster <- function(chain, threshold, initial_width = 100, dist_func = frob_norm) {
  #initial change points
  number_of_observations <- length(chain)
  cut_points <- seq(0, number_of_observations, initial_width)
  merge_bool <- T

  states <- as.character(unique(chain))

  # fit matrices
  mats <- fit_multiple_matrices(chain, change_points = cut_points)
  mats <- lapply(mats, function(x) standardize_matrix(cols = states, rows = states, mat=x))

  # merge blocks
  distance <- c()
  merge_bool <- c()
  for(i in seq(1, length(mats)-1, 2)) {
    dist <- dist_func(mats[[i]], mats[[i+1]])
    distance <- c(distance,dist)
    merge_bool <- c(merge_bool, dist < threshold)
  }
  cut_points <- cut_points[!merge_bool]
  mats <- fit_multiple_matrices(chain, change_points = cut_points)
  mats <- lapply(mats, function(x) standardize_matrix(cols = states, rows = states, mat=x))

  return( list('matrices' = mats,
               'change_points' = cut_points))
}






