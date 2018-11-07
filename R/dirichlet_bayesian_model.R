

## Priors ----

### Transition matrix of observables prior ----
obs_tmatrix_prior_sampler <- function(llambda_mat) {
  mat <- DirichletReg::rdirichlet(nrow(llambda_mat), alpha =  llambda_mat)
  colnames(mat) <- colnames(llambda_mat)
  rownames(mat) <- rownames(llambda_mat)
  mat
}


obs_tmatrix_prior_density <- function(observations, llambda_mat, log_density=F) {
  if (!log_density){
    prod(DirichletReg::ddirichlet(x = observations, alpha = llambda_mat, log = log_density))
  } else {
    sum(DirichletReg::ddirichlet(x = observations, alpha = llambda_mat, log = log_density))
  }
}


### Transition of hidden state prior
latent_state_tmatrix_sampler <- function(a,b, number_of_transitions){
  # sample betas jointly
  beta_samples <- rbeta(number_of_transitions-1, a, b)
  counterpart_beta_samples <- 1 - beta_samples

  # set diagonal entries, last entry being 1
  ret_mat <- diag(c(beta_samples, 1))

  #  modify  diagonal + 1 col
  offset_index <- seq_len(number_of_transitions-1)
  diag(ret_mat[offset_index,offset_index+ 1]) <- counterpart_beta_samples

  # return matrix
  ret_mat
}

latent_state_tmatrix_density <- function(x, a,b, log_density=F){
  if (! log_density) {
    prod(dbeta(daig(x), a, b))
  } else {
    sum(dbeta(daig(x), a, b, log = T))
  }
}

### Number of change points prior ----
number_of_cps_sampler <- function(max_num = 10) {
  sort(sample(x = seq_len(max_num), size = 1))
}

number_of_cps_density <- functionx(x, max_num = 10, log_density =F) {
  if (! log_density) {
    1/max_num
  } else {
    log(1/ max_num)
  }
}

### Sample change points given number of changes
cond_change_points_sampler <- function(number_of_change_points, index) {

  unique(c(0, sort(sample(index, size = number_of_change_points, replace = F)), length(index)))
}

cond_change_points_density <- function(number_of_change_points, index, log_density = F) {
  n <- length(index)
  k <- number_of_change_points

  if (! log_density) {
    1/ choose(n, k)
  } else {
    log(1/ choose(n, k))
  }
}

### Sample change points given max
change_point_sampler <- function(max_num, index) {
  number_of_change_points <- number_of_cps_sampler(max_num)
  cond_change_points_sampler(number_of_change_points, index)
}

change_point_density <- function(x, max_num, index, log_density = F) {
  number_of_change_points <- length(x)
  if (! log_density) {
    number_of_cps_density(number_of_change_points, log_density = log_density, max_num = max_num) *
      cond_change_points_density(number_of_change_points, index = index, log_density = log_density)
  } else {
    number_of_cps_density(number_of_change_points, log_density = log_density, max_num = max_num) +
      cond_change_points_density(number_of_change_points, index = index, log_density = log_density)
  }
}



## Posterior likelihoods
update_prior_matrix <- function(chain, prior_matrix) {
  counts <- markovchain::createSequenceMatrix(chain)
  counts <- standardize_matrix(prior_matrix,
                               rows = rownames(prior_matrix),
                               cols = colnames(prior_matrix))
  prior_matrix + counts
}

posterior_liklihood <- function(chain, transition_matrix, prior_matrix, log=F) {
  if (! log) {
    chain_likelihood(chain,
                     transition_matrix = transition_matrix,
                     log = log) *
    obs_tmatrix_prior_density(transition_matrix,
                              llambda_mat = prior_matrix,
                              log_density = log)
  } else {
    chain_likelihood(chain,
                     transition_matrix = transition_matrix,
                     log = log) +
      obs_tmatrix_prior_density(transition_matrix,
                                llambda_mat = prior_matrix,
                                log_density = log)
  }
}


posterior_obs_tmatrix_sampler <- function(chain, indicator_labels, prior_matrix) {
  # group chain intoblocks
  chain_pieces <- lapply(unique(indicator_labels), function(x) chain[indicator_labels==x])

  new_mats <- lapply(chain_pieces, function(x) update_prior_matrix(chain = x, prior_matrix = prior_matrix))
  lapply(new_mats, function(x) obs_tmatrix_prior_sampler(llambda_mat = x))
}

state_sampler <- function(chain, number_of_transitions, obs_trans_mats, state_matrix) {
  chain_index <- seq_along(chain)
  prob_matrix <- calculate_state_probs(chain, number_of_transitions, state_matrix, obs_trans_mats)

  last_state <- tail(chain_index,1)
  states <- numeric()
  states[last_state] <- number_of_transitions

  # reverse order
  for (i in last_state:2) {
    proposals <- c(states[i]-1, states[i])
    proposals <- proposals[proposals!=0]
    probs <- state_matrix[proposals, states[i]] * prob_matrix[i-1, proposals]
    probs <- probs / sum(probs)
    states[i-1] <- sample(proposals, size = 1, prob = probs)
  }

  states
}



calculate_state_probs <- function(chain, number_of_transitions, state_matrix, obs_trans_mats) {
  chain_length <- length(chain)
  transition_index <- seq_len(number_of_transitions)

  # placeholder
  indicator_probs <- matrix(NA, nrow = chain_length, ncol = number_of_transitions)

  # set initial
  indicator_probs[1,] <- 0
  indicator_probs[1,1] <- 1

  # loop forward
  for (i in 2:chain_length) {
    indicator_probs[i,] <- sapply(transition_index, function(x) {
      (t(indicator_probs[i-1,]) %*% state_matrix)[x] *
        (obs_trans_mats[[x]][chain[i-1],chain[i]])
    })
    # standardize
    indicator_probs[i,]  <- indicator_probs[i,] / sum(indicator_probs[i,])
  }
  indicator_probs
}

get_change_points <- function(states_indicators) {
  chain_index <- seq_along(states_indicators)
  change_points <- sapply(unique(states_indicators), function(x) min(chain_index[states_indicators== x]))
  change_points
}

posterior_state_sampler <- function(states_indicators, prior_a, prior_b, number_of_transitions){
  # update all but last
  freq <- table(states_indicators)
  freq <- freq[1:length(freq)-1]

  # update params
  posterior_a <- prior_a + freq
  posterior_b <- prior_b + 1

  latent_state_tmatrix_sampler(a = posterior_a, b = posterior_b, number_of_transitions = number_of_transitions)
}
