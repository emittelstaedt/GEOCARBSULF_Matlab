%% GEOCARBSULF - Set Values to Initial Condition at Simulation Start Time 
%   Filename: GEOCARBSULF_setinitialcondition.m
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
% This function sets the calculation value to the initial conditions set in
% the input files read in by GEOCARBSULF and, if desired, resampled according
% to either a Normal or Log Normal probability disribution. See other 
% functions for more details on the MC error estimation resampling.
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
%       Fluxes ("f" prefix) are in units of 10^18 mol Myrs-1
%       Rates ("k" prefix) are in units of Myrs-1
%       Stable isotopic compositions ("d" prefix) are in per mil units
%
%% Assign values  (using the loop variable irs to get proper resample run)

% atmospheric CO2 (ratio between CO2 at time t to the Pleistocene mean 
% [taken as 250 ppm]); this initial value is a place-holder 
% (it is solved for explicitly)
RCO2 = 10;

% these variables are recalculated at each time-step
oxygen = oxygen_start(irs);
Gy = Gy_start(irs);
Cy = Cy_start(irs);
Ca = Ca_start(irs);
Ssy = Ssy_start(irs);
Spy = Spy_start(irs);
dlsy = dlsy_start(irs);
dlcy = dlcy_start(irs);
dlpy = dlpy_start(irs);
dlpa = dlpa_start(irs);
dlgy = dlgy_start(irs);
dlga = dlga_start(irs);

% d34S value of old CaSO4 sulfur
dlsa = (dlst(irs)*ST(irs)-(dlpy*Spy+dlsy*Ssy+dlpa*Spa(irs)))/Ssa(irs);

% d13C value of old crustal carbonate carbon
dlca = (dlct(irs)*CT(irs)-(dlgy*Gy+dlcy*Cy+dlga*Ga(irs)))/Ca;
Rcy = Rcy_start(irs);
Rca = Rca_start(irs);

