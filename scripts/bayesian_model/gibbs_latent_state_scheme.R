

# Gibbs sampling on latent state, Bayesian model

## Initialize Model
# data
chain <- jtils::unserialize_robject('./test_data/data/test_chain.ser')

# prior param
states <- sort(unique(chain))
llambda_mat <- as.matrix(replicate(4, rep(0.25,4)), dimnames = list(states, states))
number_of_transitions <- 11
chain_index <- seq_along(chain)
transition_index <- seq_len(number_of_transitions)
chain_length <- tail(chain_index,1)
colnames(llambda_mat) <- states
rownames(llambda_mat) <- states

transition_a <- 0.89
transition_b <- 0.1

# Learned variables

state_matrix <- latent_state_tmatrix_sampler(a = transition_a, b = transition_b, number_of_transitions = number_of_transitions)
while (1 %in% diag(state_matrix)[-length(state_matrix)]) {
  state_matrix <- latent_state_tmatrix_sampler(a = transition_a, b = transition_b, number_of_transitions = number_of_transitions)
}

obs_trans_mats <- lapply(transition_index, function(x) obs_tmatrix_prior_sampler(llambda_mat))
states_indicators <- sort(sample(transition_index, replace = T, size = length(chain_index)))
get_change_points(states_indicators)

obs_trans_res <- list()
states_indicators_res <- list()
states_mat_res <- list()

# update obs matrices
for (i in 1:1000) {
  print(i)

  obs_trans_mats <- posterior_obs_tmatrix_sampler(chain, indicator_labels = states_indicators, prior_matrix = llambda_mat)
  states_indicators <- state_sampler(chain, number_of_transitions, obs_trans_mats, state_matrix)
  state_matrix <- posterior_state_sampler(states_indicators,
                                          prior_a = transition_a,
                                          prior_b = transition_b,
                                          number_of_transitions = number_of_transitions)
  obs_trans_res[[i]] <- obs_trans_mat
  states_indicators_res[[i]] <- states_indicators
  states_mat_res[[i]] <- state_matrix

  if (i %% 10 == 0) {
    plot_chain(chain, get_change_points(states_indicators))
  }

}

plot_chain(chain, get_change_points(states_indicators))


get_change_points(states_indicators_res[[100]])

chain_indicators <- jtils::unserialize_robject('./bayesian_model/states_ind.ser')
length(chain_indicators)

change_points <- lapply(chain_indicators, get_change_points)

dist <- function(x) {
  x <- x[-length(x)]
  sqrt(sum(x-c(1,1:9*1000))^2)
}

dists <- sapply(change_points, dist)
change_points[[934]]
plot_chain(chain, change_points[[934]])
library(magrittr)
