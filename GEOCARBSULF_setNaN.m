%% GEOCARBSULF - Set all Array and Constant Values for this run/time = NaN
%   Filename: GEOCARBSULF_setNaN.m
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
% Note from Royer, 2014 version of GEOCARBSULF: 
%   This function calls GEOCARBSULF_iterativeCO2solve.m to calculate 
% atmospheric CO2 (ppm) through iterative convergence (see FORTRAN scripts 
% of Park & Royer 2011 for details); this is done by calculating the 
% climatic factors that affect silicate weathering rates normalized to the 
% present-day (fwsi_climate, calculated below; called "Fbbs" in BASIC 
% scripts); fwsi_climate combines the expressions "fBt(CO2)" and "fBb(CO2)" 
% in Berner (2004) (see his equations 2.6-7, 2.29, & 5.2); through 
% inversion of fwsi_climate ("W", "V", & "X" below) and comparison to 
% fwsi_no_climate (calculated above), RCO2 can be determined; see pp. 72-76 
% in Berner (2004) for details. 
% 
%   Iterative convergence is necessary because RCO2--taken from the previous 
% time-step--is needed to calculate some of the dependent parameters; the 
% initial calculation of the new time-step for RCO2 is therefore not likely 
% to be correct
%
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
%% Set NaN values

% fill all time arrays with NAN
Sr(tt,irs) = NaN;
d13C(tt,irs) = NaN;
d34S(tt,irs) = NaN;
fR(tt,irs) = NaN;
fL(tt,irs) = NaN;
fA(tt,irs) = NaN;
fAw_fA(tt,irs) = NaN;
fD(tt,irs) = NaN;
GEOG(tt,irs) = NaN;
RT(tt,irs) = NaN;
fSR(tt,irs) = NaN;
fC(tt,irs) = NaN;

% fill all constants with NaN
ACT(irs) = NaN;
ACTcarb(irs) = NaN;
VNV(irs) = NaN;
NV(irs) = NaN;
exp_NV(irs) = NaN;
LIFE(irs) = NaN;
GYM(irs) = NaN;
FERT(irs) = NaN;
exp_fnBb(irs) = NaN;
deltaT2X(irs) = NaN;
GLAC(irs) = NaN;
J(irs) = NaN;
n(irs) = NaN;
Ws(irs) = NaN;
exp_fD(irs) = NaN;
Fwpa_0(irs) = NaN;
Fwsa_0(irs) = NaN;
Fwga_0(irs) = NaN;
Fwca_0(irs) = NaN;
Fmg_0(irs) = NaN;
Fmc_0(irs) = NaN;
Fmp_0(irs) = NaN;
Fms_0(irs) = NaN;
Fwsi_0(irs) = NaN;
Xvolc_0(irs) = NaN;
CAPd13C_0(irs) = NaN;
CAPd34S_0(irs) = NaN;
oxygen_start(irs) = NaN;
Gy_start(irs) = NaN;
Cy_start(irs) = NaN;
Ca_start(irs) = NaN;
Ssy_start(irs) = NaN;
Spy_start(irs) = NaN;
dlsy_start(irs) = NaN;
dlcy_start(irs) = NaN;
dlpy_start(irs) = NaN;
dlpa_start(irs) = NaN;
dlgy_start(irs) = NaN;
dlga_start(irs) = NaN;
Rcy_start(irs) = NaN;
Rca_start(irs) = NaN;
Rv_start(irs) = NaN;
Rg_start(irs) = NaN;
Fob(irs) = NaN;
COC(irs) = NaN;
Ga(irs) = NaN;
Ssa(irs) = NaN;
Spa(irs) = NaN;
ST(irs) = NaN;
dlst(irs) = NaN;
CT(irs) = NaN;
dlct(irs) = NaN;
kwpy(irs) = NaN;
kwsy(irs) = NaN;
kwgy(irs) = NaN;
kwcy(irs) = NaN;

