%% GEOCARBSULF - Resample All Arrays and Constants used in GEOCARBSULF 
%   Filename: GEOCARBSULF_performresampling.m
%
%   Based upon the R script written by Dana Royer for the publication Royer 
% et al (2014). Please see the ROYER_README.txt file in the same directory
% as this script for details on their work and the R code. 
%
%   This script translates a portion of the R code to Matlab with a
% reorganization into functions to make the overall code structure more
% readable and easier to understand than the single file R code.
% 
% See the GEOCARBSULF_main.m comments for more detailed notes
%
% This function calls GEOCARBSULF_resampletimearray.m to resample the 
% time arrays and GEOCARBSULF_resampleconstants.m to resample constant 
% values across resampleN model calculations.  Note constants are still 
% held constant for each model simulaton, but vary according to their 
% resampled distribution between simulations. This allows for the Monte 
% Carlo based method of esimating uncertainty in the solutions.
%
% Note: this is not a function as it requires passing in and out too much
% information to the main code. 
%
% Some additional notes on variables and abbreviations used
%       y=young;            a=old; 
%       p=pyrite;           s=sulfate; 
%       c=carbonate;        si=silicates; 
%       g=organic matter;   b=burial; 
%       m=degassing;        w=weathering
%
%   UNITS
%       Masses are in units of 10^18 mol
%       Fluxes ("F" prefix) are in units of 10^18 mol Myrs-1
%       Rates ("k" prefix) are in units of Myrs-1
%       Stable isotopic compositions ("d" prefix) are in per mil units
%
% ***************************  TIME ARRAYS ********************************












% ****************************  CONSTANTS  ********************************












