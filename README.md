# Folder structure

## Master_Thesis_ZAETTA_Paul.pdf

This file contains the master thesis, which is "Very short-term analysis of wind power generation in a probabilistic forecasting framework". 

## Data_2016.xlsx

It contains the final data set. It has been previously cleaned using python (see the python script). The data set for wind power
generation in Galicia consists of 52123 valid data points (from 1th January 2016 to 31th December 2016).

## Cross_Validation_code_Matlab

It contains all the code (Matlab) necessary to use the cross validation method for adjusting and calculating predictions to estimate the performance of the prediction model. More precisely, it helps us to determine the optimal parameters of the models. 

## MAE_RMSE_code_Matlab

It contains all the code (Matlab) necessary to use the Mean Average Error and Root Mean Square Error methods. They allow us to evaluate the point forecasts by measuring differences between observed and predicted values.

Warning: the code "MAE_RMSE_ARMA_and_final_chart" contains both the MAE and RMSE levels and also the final figure, which represents a prediction sequence of wind generation using ARMA dynamics (page 27 in the "Master_Thesis_ZAETTA_Paul.pdf").

## CRPS_code_Matlab

It contains all the code (Matlab) necessary to use the Continuous Rank Probability Score. It is a verification tool related to probabilistic forecasting systems, and it is a quantity that highlights both the forecast of sharpness and calibration.

## functions_code_Matlab

It contains the other functions that were used during this study. The GL_transform function corresponds to the generalised logit transform and the IGL_transform corresponds to the invers generalised logit transform. 
