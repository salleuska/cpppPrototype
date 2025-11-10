## checks
library(cpppPrototype)
x <- new_cppresults(cppp = 0.4, ppp = runif(10), obs_ppp = 0.37)
class(x)          # should show "cpppResult" "list"
length(x$ppp)     # should be 10
