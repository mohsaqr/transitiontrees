# Cap threads (CRAN: CPU time close to elapsed time).
Sys.setenv(OMP_THREAD_LIMIT = 2)

library(testthat)
library(transitiontrees)

test_check("transitiontrees")
