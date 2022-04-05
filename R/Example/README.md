Details on key files in the folder 'Example':

**data_gen.R**

Generates VAR transition matrix using function 'Gen_A' and generates true error covariance matrix using function 'Gen_Sigma' for the model by Ghysels as in the paper 'Macroeconomics and the reality of mixed frequency data' and then generates data from the model.

**BMF.R**

This code first generates data from the model by Ghysels using 'data_gen.R' for setting 1 in the paper, fit the proposed BMF model to the simulated data and then provides RMSE, CRPS and log-score values.

**MFBVAR.R**

Code for implementing the estimation by R package `mfbvar', using the model by Schorfheide and Song as in the paper 'Real-Time Forecasting With a Mixed-Frequency VAR' by Schorfheide and Song.

**midas.R**

Code for fitting MIDAS regression model by R package `midasr', as given in the paper 'Midas regressions: Further results and new directions' by Ghysels, Sinko and Valkanov.

**benchmark.R**

Code for fitting the benchmark model (a random walk model with drift) used in forecasting for obtaining the relative RMSE values.

**functions.R**

Contains all the R functions required to fit the proposed BMF model as well as other competing models.
