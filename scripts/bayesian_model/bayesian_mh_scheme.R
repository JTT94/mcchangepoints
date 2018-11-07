nsim <- 10^5
prior_mat <- attr(markovchainFit(chain)$estimate,'transitionMatrix')
sample_model <- function(chain, prior_mat, number_of_points) {
  chain_index <- seq_along(chain)
  # sample params

  # number of transitions

  num_changes <- number_of_points  + sample(c(-1,0,1),size = 1)
  num_changes <- max(0, num_changes)

  # change points
  change_points <- cond_change_points_sampler(index = chain_index,
                                              number_of_change_points = num_changes)

  # sample matrices
  tran_mats <- lapply(seq_len(num_changes+1), function(x) obs_tmatrix_prior_sampler(llambda_mat = prior_mat))

  # cut index
  cut_labels <- cut(chain_index, breaks = change_points, labels = F)
  chain_pieces <- lapply(unique(cut_labels), function(x) chain[cut_labels==x])

  # calculate likelihood
  likelihood <- sum(sapply(seq_along(chain_pieces), function(x) posterior_liklihood(chain = chain_pieces[[x]],
                                                                                   transition_matrix = tran_mats[[x]],
                                                                                   prior_matrix = prior_mat,
                                                                                   log = T)))

  return(list(transition_matrices = tran_mats,
              change_points = change_points,
              likelihood = likelihood))
}

last_state<- sample_model(chain, prior_mat)
best_fit <- initial_estimates
for ( i in seq_len(nsim)) {
  print(i)
  prop <- sample_model(chain, prior_mat, number_of_points = length(last_state$change_points))

  if (prop$likelihood > best_fit$likelihood) {
    best_fit <- prop
  }

  if (log(runif(1)) < (prop$likelihood - last_state$likelihood)){
    last_state <- prop
  }
}


par(mfrow=c(1,1))
plot_chain(chain, change_points = best_fit$change_points)
