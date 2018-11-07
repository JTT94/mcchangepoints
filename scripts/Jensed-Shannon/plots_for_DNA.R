dat <- read.fasta(file = "../changepoints/data/chrI.fa.gz")
dat2 <- unlist(dat, use.names = FALSE)
dat2f <- as.factor(dat2)[1:1500]

DNA_plot(dat2f, c(300))

perform_DNA_test(JSD_algorithm, size = 1500, cutoff=15)

cps = find_change_points_JSD(dat2f)

DNA_fractions(dat2f, cps[-1])
