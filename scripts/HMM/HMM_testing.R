library(seqinr)
library(zoo)
library(ggplot2)
library(dplyr)

#==============================================================================
# Testing on data
#==============================================================================

#==================================================================
# Load simulated data

d <- as.numeric(jtils::unserialize_robject('../test_data/data/test_chain.ser'))[1:10000]
obs <- sort(unique(d))
n_obs <- length(obs)

#------------------------------------------------------------------------------

# Load real data
dat <- read.fasta(../data/file="../data/chrI.fa.gz")
dat2 <- unlist(dat, use.names=FALSE)
dat2[dat2=="a"] <- 1
dat2[dat2=="c"] <- 2
dat2[dat2=="g"] <- 3
dat2[dat2=="t"] <- 4
dat2 <- as.numeric(dat2)

d <- dat2[1:10000]

#------------------------------------------------------------------------------

# Construct initial matrix A for hidden state transitions
# Give initial guess of max K, and lambda
cpt_number <-10
lambda <- 0.001

A <- matrix(0, nrow=cpt_number, ncol=cpt_number)
A[row(A)-col(A) == -1] <- lambda
A[row(A)-col(A) == 0] <- 1-lambda
A[row(A) + col(A) == 2*cpt_number] <- 1

# Construct pi
pi <- c(1, rep(0, cpt_number-1))

# Construct initial matrix B for emmission probabilities
B <- array(1/n_obs, c(n_obs, n_obs, cpt_number))


#------------------------------------------------------------------------------

# Carry out the parameter fitting
bw <- BW_recursion(d, pi, A, B, max_iter=20)
A_new <- bw$A_update
B_new <- bw$B_update
pi_new <- bw$pi_update

#------------------------------------------------------------------------------

# Calculate most likely state path
path <- viterbi(d, pi_new, A_new, B_new)


# Produce state path plot similar to Figures 1 and 2
d2 <- d
d2[d2==1] <- "A"
d2[d2==2] <- "C"
d2[d2==3] <- "G"
d2[d2==4] <- "T"

path_plot <- as.data.frame(d2) %>%
  bind_cols(as.data.frame(path)) %>%
  mutate(n = 1:n()) %>%
  rename(Nucleotide = d2)

plot1 <- ggplot(as.data.frame(path_plot), aes(n, path))
plot1 +
  geom_point(aes(colour = Nucleotide), shape=124, size=6) +
  theme_bw() +
  scale_y_continuous(name="Segment",breaks=seq(1,cpt_number,1), minor_breaks = NULL) +
  scale_x_continuous(name="Position",breaks=seq(0,length(d),1000), minor_breaks = NULL)

#------------------------------------------------------------------------------

# Produce plot of posterior probability of being in a particular state at any time
f <- log_forward_alg(d, pi_new, A_new, B_new)
b <- log_backward_alg(d, A_new, B_new)
logdenom <- log_denom(f)
gam <- log_gamma(f,b, logdenom)
plot(exp(gam)[6,], type='l')

