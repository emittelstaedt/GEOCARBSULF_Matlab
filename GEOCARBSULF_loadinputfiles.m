%% GEOCARBSULF - Read Inputs and Prepare Output Files
%   Based upon the R script written by Dana Royer for the publication Royer 
% et al (2014). Please see the ROYER_README.txt file in the same directory
% as this script for details on their work and the R code. 
%
%   This script translates a portion of the R code to Matlab with a
% reorganization into functions to make the overall code structure more
% readable and easier to understand than the single file R code.
% Additionally, a time-varying degassing term Xvolc was added to allow
% calculation of the CO2 contribution of hotspot volcanic fluxes
%
% Note that since the input file names do not change, they are hard coded
% here and assumed to reside within the same directory as this code. 
%
%% Load Input Files

%reading in the two input files (see ROYER_README.txt and noted papers)
inputs = readtable("GEOCARB_input_summaries.csv");
time_arrays = readtable("GEOCARB_input_arrays.csv");

% determine number of time steps
nsteps = height(time_arrays);

% Check if Monte Carlo or Individual Resample
if (strcmp(loop_parameters,'TRUE'))
    % resample one parameter at a time
    % when start loop, for main code, one parameter will be set to "TRUE"
    % for each resample loop.
    inputs(:,3) = 'FALSE';
else
    % running monte carlo resample routines
    % if resampleN set = 1, will not resample at all
end










