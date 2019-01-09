# MGN
Modified Gauss-Newton algorithm for Hankel Structured Low-rank Approximation

# The R code used in paper 

**Please install auxiliary binary package from [rhlra2](rhlra2/) subdirectory before running examples in R**.

Following public packages are also needed for installing:
`install.packages(c("fftw", "svd", "Matrix", "orthopolynom"))`

## Directories description

### [auxiliary](auxiliary/)

Contains an implementation of used method (MGN and VPGN)

### [plots](plots/)

Contains code for creating plots printed in paper (comparison of MGN/VPGN methods by
accuracy and speed).

### [examples](examples/)

Contains some simple examples of methods usage.

## Common rule for all R files

You need to set a working directory (either in RStudio or using `setwd`) to source file,
then run all lines from it. Needed dependencies will be evaluated automatically.

To begin with, you can run examples from [examples](examples/) directory. It contains description
on how to run MGN/VPGN method.