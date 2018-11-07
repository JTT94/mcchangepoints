# Read genomic data

dat <- read.fasta(file = "../changepoints/data/chrI.fa.gz")
dat2 <- unlist(dat, use.names = FALSE)
dat2f <- as.factor(dat2)

serial_time = numeric(10)
parallel_time = numeric(10)

# simulate for different sizes to understand time complexity

for (size in 6:10)

{
  print(size)

  test = dat2f[200:(3000*size)]

  result_serial = system.time(results <- JSD_algorithm(test, cutoff = 15))
  serial_time[size] = result_serial[3]

  print(result_serial)


  result_parallel = system.time(results <- JSD_algorithm(test, cutoff = 15, par = TRUE))
  parallel_time[size] = result_parallel[3]

  print(result_parallel)
}

length_list = 3000*(1:10) - 200

plot(length_list, parallel_time, col="blue", type="line", ylab="time", xlab="chain size")
lines(length_list, serial_time, col="red")
legend("topleft", col=c("blue", "red"), legend=c("parallel", "serial"), lty=1)
