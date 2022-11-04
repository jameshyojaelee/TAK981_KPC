library(doParallel)
library(foreach)

starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(1)  # use 1 core, set to the number of our cores
registerDoParallel(4)  # use 1 core, set to the number of our cores
registerDoParallel(numCores)  # use multicore, set to the number of our cores
