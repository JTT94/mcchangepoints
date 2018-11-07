# Markov chain helpers

fit_transition_matrix <- function(chain) {
  attr(markovchainFit(chain)$estimate,'transitionMatrix')
}

count_transitions <- function(chain){
  tab <- table(sapply(seq_len(length(chain)-1), function(x) paste0(chain[x],"->",chain[x+1])))
  tab
}

# helper to ensure matrix is complete in the case where 0 transitions occur
standardize_matrix <- function(rows, cols, mat){

  mat <- as.data.frame(mat)
  for (row in rows) {
    if (! row %in% rownames(mat)) {
      mat[row,] <- 0
    }
  }
  for (col in cols) {
    if (! col %in% colnames(mat)) {
      mat[,col] <- 0
    }
  }
  as.matrix(mat)
}


# distance functions
frob_norm <- function(A,B) {
  sqrt(sum((A-B)^2))
}


#likelihood of sequence given transition_matrix
chain_likelihood <- function(chain, transition_matrix, log = F) {
  update_freq <- createSequenceMatrix(chain)
  update_freq <- standardize_matrix(update_freq,
                                    cols = colnames(transition_matrix),
                                    rows = rownames(transition_matrix))
  if (! log) {
    prod(transition_matrix ^ update_freq)
  } else {
    sum(log(transition_matrix) * update_freq)
  }

}

#plot cumulative totals
plot_chain <- function(chain, change_points = 0) {
  # running total frequency
  totals <- data.frame(index = seq_along(chain), chain = chain, value = 1) %>%
    tidyr::spread(key = 'chain', value = 'value', fill=0) %>%
    dplyr::select(-index) %>%
    dplyr::mutate_all(.funs=cumsum)

  # sense check plot
  x <- seq_along(totals[,1])
  y <- totals
  matplot(x, y, type='l')
  abline(v=change_points)
}
