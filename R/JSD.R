library(markovchain)
library(DescTools)

# These functions implement change point detection using Markov Models for segmentation based on the Jensed-Shannon distance,
# described in the paper:
#
# Thakur, Vivek, Rajeev K. Azad, and Ram Ramaswamy. "Markov models of genome segmentation." Physical Review E 75.1 (2007): 011915.
#
# Up to now, only the model selection framework is implemented for the statistical significance test (Sec. II.A of the previous).
#
# The main function is the JSD_algorithm one


#' Entropy of a MC chain
#'
#' @description Computes the entropy of a sequence, assuming a MC structure.
#'
#' @param  sequence_matrix = the observed sequence matrix of which to estimate the entropy.
#'
#' @return  The computed value of entropy
#'
#' @export
#'

compute_entropy_MC = function(sequence_matrix) {
  # function computing entropy for a MC sequence.

  occurrences = rowSums(sequence_matrix)  # occurrences for each letter

  transition_probabilities = sequence_matrix / replace(occurrences, occurrences ==
                                                         0, 1)
  # transition_probabilities[i,j] contains the prob of going to j starting from i. Replace is used to avoid 0/0

  letter_probabilities = occurrences / sum(sequence_matrix)  # in this formulation, this is the same as (length(chain) - 1)

  entropy = -sum(letter_probabilities * rowSums(transition_probabilities *
                                                  replace(
                                                    log2(transition_probabilities),
                                                    log2(transition_probabilities) == -Inf,
                                                    0
                                                  )))

  # replace is to avoid NA values that can be caused by 0*Inf. Remember that in computation of entropy 0*log2(0)=0

  return (entropy)
}


#' Find best split of a sequence
#'
#' @description Find the best possible split of a sequence according to the Jensen-Shannon divergence. It computes this measure for all possible split points
#' and select the best split as the one maximizing it
#'
#' @param  chain = the sequence of which to find the split point
#' @param  cutoff = smaller possible size of a chunk

#' @return  a list containing i) index of the change point ii) value for the Jensen-Shannon divergence at that point (used in the BIC criterion).
#' @export

find_best_split = function(chain, cutoff = 15) {
  # we find the best split for the chain as the split_point maximing the D_max. Could be parallelized?

  # find the states of the MC in the followin way
  states = rownames(createSequenceMatrix(chain))

  # now define the starting matrices:
  sequence_matrix_1 = createSequenceMatrix(chain[1:(cutoff-1)], possibleStates = states)
  sequence_matrix_2 = createSequenceMatrix(chain[(cutoff-1) : length(chain)], possibleStates = states)


  #' JS divergence of two subsequences. It is defined inside `find_best_split` so that the variables `sequence_matrix_1` and `sequence_matrix_2`
  #' are visible to it.
  #'
  #' @description Computes the JS divergence of the two adjacent subsequences delimited by `split_point`. Actually, the entroy of the overall sequence should be added
  #' as well, but is neglected since it is the same for all the possible split points.
  #'
  #' @param  chain = the sequence of which to estimate the JS divergence
  #' @param  split_point = the position in which the sequence is split
  #' @param  w = 2-elements vector containing the weights to apply to the two different subsequences in the computation. If `NULL`, the weights are taked to be
  #' proportional to the length of the subsequence

  #' @return  The computed value of JS divergence
  #' @export

  compute_JS_div = function(chain, split_point, w = NULL) {

    if (is.null(w))
      w = c(split_point / length(chain), 1 - split_point / length(chain))
    else
      w = w / sum(w)  # normalize the weights

    sequence_matrix_1[chain[split_point - 1], chain[split_point]] <<- sequence_matrix_1[chain[split_point - 1], chain[split_point]] + 1
    sequence_matrix_2[chain[split_point - 1], chain[split_point]] <<- sequence_matrix_2[chain[split_point - 1], chain[split_point]] - 1

    #print(compute_entropy_MC(sequence_matrix_1))

    return(-w[1] * compute_entropy_MC(sequence_matrix_1) - w[2] * compute_entropy_MC(sequence_matrix_2))
  }


  D_list = sapply(seq(cutoff, length(chain) - cutoff), function(index, chain)
    compute_JS_div(chain, index), chain = chain)

  D_max = max(D_list)
  best_index = seq(cutoff, length(chain) - cutoff)[which.max(D_list)]

  return(list("index" = best_index, "D" = D_max))
}



#' Jensen-Shannon divergence algorithm
#'
#' @description Recursive function. In a top down manner, it finds one change point for the whole sequence and then performs the same operation
#' on the two resulting subsequences. It stops according to a model selection criterion, or when the smallest size of the sequence is reached.
#'
#' @param  chain = the (supposedly) piecewise Markov Chain on which to run the algorithm
#' @param  cutoff = smaller possible size of a chunk
#' @param  t_in = the starting position of the considered subsequence in the current call of the function. In the first call, `t_in = 1`` if you want to
#' perform the algorithm on the whole sequence
#' @param  t_fin = the final position of the considered subsequence in the current call of the function. In the first call, `t_fin = NULL` if you want to
#' perform the algorithm on the whole sequence. This forces t_fin = length(chain) in the function
#'
#' @return  a 1-row matrix containing vector of change points. The first element of the matrix is always `1``, due to the recursive structure
#' @export

find_change_points_JSD = function(chain,
                                  cutoff = 15,
                                  t_in = 1,
                                  t_fin = NULL) {
  if (is.null(t_fin))
    t_fin = length(chain)

  if (t_fin - t_in < 2 * cutoff)
    return(t_in)


  res = find_best_split(chain[seq(t_in, t_fin)], cutoff = cutoff)

  new_split = res$index + t_in
  D = res$D + compute_entropy_MC(createSequenceMatrix(chain[seq(t_in, t_fin)]))  # JS distance. Add entropy of whole chain as the inside
  # function only computes the one of subsequences

  # the following significance condition comes from the BIC criterion
  N = t_fin - t_in + 1
  if (2 * D * N < 16 * log(N))
    # the hard-coded factor 16 comes from assuming a 1-th order Markov process
    return(t_in)

  return (cbind(
    Recall(
      chain,
      cutoff = cutoff,
      t_in = t_in,
      t_fin = new_split
    ),
    Recall(
      chain,
      cutoff = cutoff,
      t_in = new_split,
      t_fin = t_fin
    )
  ))
}



#' Jensen-Shannon divergence based Markov Model for segmentation
#'
#' @description This is a wrapper for the algorithm. After having computed the change points, it estimates the matrices in the varios parts.
#'
#'
#'
#' @param  chain = the (supposedly) piecewise Markov Chain on which to run the algorithm
#' @param  cutoff = smaller possible size of a chunk
#'
#' @return  a list containing i) vector of change points ii) estimated transition matrices for each chunk
#' @export

JSD_algorithm = function(chain, cutoff) {

  change_points = find_change_points_JSD(chain, cutoff = cutoff)

  change_points = change_points[1, seq(2, length(drop(change_points)))]  # discard useless dimension and the first element, that is always =1. This may give error if no change point is found

  mle_mats = fit_multiple_matrices(chain, change_points)

  return(list(
    "change_points" = change_points,
    "estimated_matrices" = mle_mats
  ))

}
