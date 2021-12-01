# Paper replication: simplified version of Khan & Thomas (2007)

This is the second homework of QMM II course. The matlab file is organised in the following way:

1. main.m is the main file to run. I did not use separate function to find price. Since direct iteration using bisection methods can handle the problem well.
2. func.m collects all miscellaneous functions that are used in main.m
3. redundant.m collects some code that is not used in the main.m but may becomes useful.
4. latex.m is the function to produce latex table from Matlab.

# folders
1. Graph folder and latex folder (not show up here) are independent of the main.m
2. Fortran folder is the Fortran code of the orginal paper.

# Current problems
Market clear tolerance can only be set at 1e-6 in order to converge. It could be due to:

a) I use eigenvector to find stationary distribution of firms. 

b) the interpolation method is spline but slightly different from what is given in problem set. 

c) I find some values are slightly negative (e.g. -0.0000001), and I force them to be zero. and sometimes firms are willing to delete inventory. 

d) I attempted to solve this problem using fzero, linear approxiation, quasi newton (some code provided in redundant.m), but the result is pretty much the same. But price_new - price_old do not have opposite sign at two bounds of the price (3.2 and 3.3). This means that we need to narrow the search range before using optimizer.


## Github link
- https://github.com/yifan-lyu/QMM_hw2

