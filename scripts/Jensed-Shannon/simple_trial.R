# load test data. You have to install the package jtils:
# devtools::install_github('JTT94/jtils')

chain <-
  jtils::unserialize_robject('./test_data/data/test_chain.ser')

cutoff = 15  # cutoff of smallest possible chunk

results = JSD_algorithm(chain, cutoff)

change_points = results$change_points
matrices = results$estimated_matrices

# we could devise some metric of closeness between matrices in order to assess the quality of the algorithm.

original_matrices = jtils::unserialize_robject('./test_data/data/test_tmatrices.ser')

norms = numeric(10)

for (i in 1:10)
{
  norms[i] = norm(matrices[[i]][, ] - original_matrices[[i]], type = "F")  # simple measure of similarity
}
summary(norms)
