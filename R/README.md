
Details on key files in the folder 'R':

**mat_gen_BMF.R**

Generates VAR transition matrix using function 'Gen_A' and generates true error covariance matrix using function 'Gen_Sigma' for the proposed BMF model.

**data_gen_neutral_1.R**

Generates VAR transition matrix using function 'Gen_A' and generates true error covariance matrix using function 'Gen_Sigma' for the model by Ghysels as in the paper 'Macroeconomics and the reality of mixed frequency data' and then generates data from the model.

**data_gen_neutral_2.R**

Generates VAR transition matrix using function 'Gen_A' and generates true error covariance matrix using function 'Gen_Sigma' for the model by Schorfheide and Song as in the paper 'Real-Time Forecasting With a Mixed-Frequency VAR' and then generates data from the model.

**forecast_BMF.R**

This code first generates data from the proposed BMF model, estimates the model using Gibbs sampling algorithm, provides forecasts and forecast errors and also provides the metrics regarding the model selection and estimation consistency of the BMF model.

**mfbvar.R**

Code for implementing the estimation by R package `mfbvar', using the model by Schorfheide and Song as in the paper 'Real-Time Forecasting With a Mixed-Frequency VAR' by Schorfheide and Song.

**midas.R**

Code for fitting MIDAS regression model by R package `midasr', as given in the paper 'Midas regressions: Further results and new
directions' by Ghysels, Sinko and Valkanov.

**random_walk.R**

Code for fitting the benchmark model used in forecasting : a random walk model with drift.

**quarterly_VAR.R**

Code for fitting a quarterly VAR model which was fitted to the real macroeconomic data used in the paper.

**prepare_data1.R**

Code to prepapre Data 1 (as mentioned in the paper), the first subset of the real macroeconomic data.

**prepare_data2.R**

Code to prepapre Data 2 (as mentioned in the paper), the second subset of the real macroeconomic data.

**prepare_data3.R**

Code to prepapre Data 3 (as mentioned in the paper), the entire real macroeconomic data.
