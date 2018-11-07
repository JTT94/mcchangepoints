
# Load requirements
library(DirichletReg)
library(markovchain)

library(tidyr)
library(magrittr)
library(dplyr)

# Genrate random transition matrix, rows sum to 1
gen_transition_matrix <- function(number_of_states){
  random_matrix <- DirichletReg::rdirichlet(number_of_states, alpha = rep(1/number_of_states, number_of_states))
  rownames(random_matrix) <- as.character(seq_len(number_of_states))
  colnames(random_matrix) <- as.character(seq_len(number_of_states))

  random_matrix
}

gen_latent_state_transition_matrix <- function(a,b, number_of_transitions){
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

gen_initial_state <- function(number_of_states){
  initial_state <- rep(0, number_of_states)
  position <- floor(10*runif(1, 0.1,0.1+number_of_states/10))
  initial_state[position] <- 1
  initial_state
}

gen_mc <- function(number_of_states, transition_matrix){
  mc_obj <- new('markovchain',
                states = as.character(seq_len(number_of_states)),
                transitionMatrix = transition_matrix
  )
}

gen_chain <- function(chain_length, mc_obj, initial_state){
  rmarkovchain(n = chain_length, object = mc_obj, t0 = initial_state)
}

# Generate Markov Chain
number_of_states <- 4
number_of_transitions <- 10
length_of_piece <- 1000

transition_matrices <- lapply(seq_len(number_of_transitions), function(x) gen_transition_matrix(number_of_states))
markov_chains <- lapply(seq_len(number_of_transitions), function(x) gen_mc(number_of_states, transition_matrices[[x]]))

initial_state <- '1'
sim_chain <- c()
for (i in seq_len(number_of_transitions)) {
  # Get MC
  mc_obj <- markov_chains[[i]]

  # Simulate chain
  sim_chain <- c(sim_chain, gen_chain(length_of_piece, mc_obj,initial_state))

  # update initial state
  initial_state <- tail(sim_chain, 1)
}

# running total frequency
totals <- data.frame(index = seq_along(sim_chain), chain = sim_chain, value = 1) %>%
  tidyr::spread(key = 'chain', value = 'value', fill=0) %>%
  dplyr::select(-index) %>%
  dplyr::mutate_all(.funs=cumsum)

# sense check plot
x <- seq_along(totals[,1])
y <- totals
par(mfrow=c(2,1))
matplot(x, y, type='l')
matplot(x, y, type='l')
abline(v=1000*1:10)

# Compare MLE
markovchain::markovchainFit(chain[1:1000])$estimate
transition_matrices[[1]]

