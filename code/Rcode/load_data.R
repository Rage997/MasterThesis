# See doc: https://cran.r-project.org/web/packages/reticulate/index.html

library(reticulate)
np <- import("numpy")
# data reading
mat <- np$load("../matrix.npy")
X <- mat

n <- dim(X)[2]
# p, d hardcoded in python
p <- 16
d <- 20
d.true <- d
