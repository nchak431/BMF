# BMF

This contains all the R codes used for implementing the methods proposed in the paper 'A Bayesian framework for sparse estimation in high-dimensional mixed frequency Vector Autoregressive models' by Chakraborty, Khare and Michailidis.

Details on the key files:

**functions.R**

Contains all the R functions related to the estimation of proposed BMF model. The main R code 'forecast_BMF.R' that uses these functions is provided in the 'R' folder.

**Rcpp_functions.cpp**

Contains the C++ functions for estimating the proposed BMF model using Rcpp which leads to significantly less computational time.

**Rcpp_main.R**

The main R code that uses 'Rcpp_functions.cpp' for estimation of the proposed BMF model using Rcpp.
