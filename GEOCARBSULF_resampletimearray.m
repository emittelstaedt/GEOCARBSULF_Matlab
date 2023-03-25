%% GEOCARBSULF - Resample Time Arrays 
%   Filename: GEOCARBSULF_resampletimearray.m
%
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
% This function generates resamples for time-dependent arrays, including 
% the potential clipping of distributions for Monte Carlo based method of
% esimating uncertainty in the solutions.
%
function [array_resamp] = GEOCARBSULF_resampletimearray(x,xsd,name,inputs,resampleN)
%
%   Function to resample the time varying arrays in the inputs files using
% either a normal or log-normal distrubution as determined by the input
% files. By default the input files have all arrays resamples using a
% gaussian normal distribution except for deltaT2X. 
 
% set up array for resample. 
array_resamp = zeros(height(x),resampleN);

% parse parameters from input summaries file
row = find(strcmp(inputs.parameter,name)==1);
resampleYorN = inputs{row,'resample'};
distrib = inputs{row,'distribution_type'};
low_lim = inputs{row,'lower_limit'};
upper_lim = inputs{row,'upper_limit'};

% if resampling=1 or if resampling is turned off for the time-dependent 
% array in question (a full matrix is generated but with resampleN mean 
% values at each time-step)
if (resampleN == 1 || strcmp(resampleYorN,'FALSE') )
    array_resamp = x{:}*ones(1,resampleN);
end

% Resample using a Normal Distribution
if (strcmp(resampleYorN,'TRUE') && (resampleN > 1))
    mean_dist = x.(name);
    std = xsd.(['e',name])./2;
    
    % Select probability distribution for resample
    if (strcmp(distrib,'gaussian'))
        temp_resample = random('Normal',mean_dist,std);

    elseif strcmp(distrib,'lognormal')
        temp_resample = random('LogNormal',log(mean_dist),log(2.*std)./2);
        
    else
        error('Selected probability distribution not included.');
    end

    % check for passing allowable limits of parameters
    temp_resample(temp_resample<=low_lim) = low_lim + 0.0001;
    temp_resample(temp_resample>=upper_lim) = upper_lim - 0.0001;

    % expand out value for all resample runs (e.g., resampleN times)
    array_resamp = temp_resample(:)*ones(1,resampleN);
end





