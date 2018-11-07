# experiment on full chromosome:

dat <- read.fasta(file = "../changepoints/data/chrI.fa.gz")
dat2 <- unlist(dat, use.names = FALSE)
dat2f <- as.factor(dat2)

time_study_full = system.time(results_full <- JSD_algorithm(dat2f, cutoff = 15))

saveRDS(list("time"=time_study_full, 'results' = results_full), file = "results_full.rds")




readRDS(file = "results_full.rds")

