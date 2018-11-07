system.time(find_best_split_parallel(chain))
system.time(find_best_split(chain))

system.time(find_change_points_JSD(chain, par=TRUE))
system.time(find_change_points_JSD(chain))
