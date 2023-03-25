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
%       Fluxes ("f" prefix) are in units of 10^18 mol Myrs-1
%       Rates ("k" prefix) are in units of Myrs-1
%       Stable isotopic compositions ("d" prefix) are in per mil units
%
% ***************************  TIME ARRAYS ********************************

% generate random distributions for the time-dependent parameters 
% alphabetical order 

% 87Sr/86Sr of shallow-marine carbonate ([87Sr/86Sr - 0.7] x 10^4)
Sr = GEOCARBSULF_resampletimearray(time_arrays(:,'Sr'),time_arrays(:,'eSr'),'Sr',inputs,resampleN);

% d13C of shallow-marine carbonate (per mil); called "DLCOC" in BASIC code
d13C = GEOCARBSULF_resampletimearray(time_arrays(:,'d13C'),time_arrays(:,'ed13C'),'d13C',inputs,resampleN);

% d34S of marine sulfate sulfur (per mil) (from Wu et al, 2010); called "DLSOC" in BASIC code
d34S = GEOCARBSULF_resampletimearray(time_arrays(:,'d34S'),time_arrays(:,'ed34S'),'d34S',inputs,resampleN);

% effect of relief on chemical weathering at time (t) relative to the 
% present-day; calculated from equation (5) in Berner (2006b)
fR = GEOCARBSULF_resampletimearray(time_arrays(:,'fR'),time_arrays(:,'efR'),'fR',inputs,resampleN);

% land area covered by carbonates at time (t) relative to the present-day
fL = GEOCARBSULF_resampletimearray(time_arrays(:,'fL'),time_arrays(:,'efL'),'fL',inputs,resampleN);

% If using Godderis et al., 2012 values for various parameters 
if (strcmp(Godderis,'TRUE'))

    % land area at time (t) relative to the present-day
    fA = GEOCARBSULF_resampletimearray(time_arrays(:,'fA_Godderis'),...
             time_arrays(:,'efA_Godderis'),'fA_Godderis',inputs,resampleN);

    % fraction of land area experiencing chemical weathering (runoff > 0)
    fAw_fA = GEOCARBSULF_resampletimearray(time_arrays(:,'fAw_fA_Godderis'),...
             time_arrays(:,'efAw_fA_Godderis'),'fAw_fA_Godderis',inputs,resampleN);

    % change in global river runoff at time (t) relative to the present-day 
    % in the absence of changes in solar luminosity and CO2 (i.e., mainly 
    % due to changes in paleogeography)
    fD = GEOCARBSULF_resampletimearray(time_arrays(:,'fD_Godderis'),...
             time_arrays(:,'efD_Godderis'),'fD_Godderis',inputs,resampleN);

    % change in land mean surface temperature for areas experiencing 
    % chemical weathering (runoff > 0) at time (t) relative to the present-
    % day in the absence of changes in solar luminosity and CO2 (i.e., 
    % mainly due to changes in paleogeography) (K)
    GEOG = GEOCARBSULF_resampletimearray(time_arrays(:,'GEOG_Godderis'),...
             time_arrays(:,'eGEOG_Godderis'),'GEOG_Godderis',inputs,resampleN); 

else % not using Godderis et al., 2012

    % VARIABLE DEFINITIONS SAME AS ABOVE
    fA = GEOCARBSULF_resampletimearray(time_arrays(:,'fA'),...
             time_arrays(:,'efA'),'fA',inputs,resampleN);

    fAw_fA = GEOCARBSULF_resampletimearray(time_arrays(:,'fAw_fA'),...
             time_arrays(:,'efAw_fA'),'fAw_fA',inputs,resampleN);

    fD = GEOCARBSULF_resampletimearray(time_arrays(:,'fD'),...
             time_arrays(:,'efD'),'fD',inputs,resampleN);

    GEOG = GEOCARBSULF_resampletimearray(time_arrays(:,'GEOG'),...
             time_arrays(:,'eGEOG'),'GEOG',inputs,resampleN);

end

% coefficient relating continental runoff to temperature change 
% (runoff/runoff(0)=1+RT*(T-T0)), as determined from the GCM simulations of 
% Godderis et al 2012 (1/K); RT is called "Y" in Berner 2004 and "RUN" in 
% GEOCARB III; in the BASIC scripts, RT is assigned a value of 0.045 during 
% times with large ice sheets and a value of 0.025 for all other times; 
% there is little difference in estimated CO2 between these two approaches 
% for RT
RT = GEOCARBSULF_resampletimearray(time_arrays(:,'RT'),time_arrays(:,'eRT'),'RT',inputs,resampleN);

% seafloor creation rate at time (t) relative to the present-day
fSR = GEOCARBSULF_resampletimearray(time_arrays(:,'fSR'),time_arrays(:,'efSR'),'fSR',inputs,resampleN);

% effect of carbonate content of subducting oceanic crust on CO2 degassing 
% rate at time (t) relative to the present-day
fC = GEOCARBSULF_resampletimearray(time_arrays(:,'fC'),time_arrays(:,'efC'),'fC',inputs,resampleN);


% ****************************  CONSTANTS  ********************************

% activation energy (E) for dissolution of Ca- and Mg-silicates on land, 
% where ACT = E/RTT0 (1/K); called "Z" in Berner (2004)
ACT = GEOCARBSULF_resampleconstant('ACT',inputs,resampleN);

% activation energy (E) for dissolution of carbonates on land (see p. 53 in 
% Berner 2004)
ACTcarb = GEOCARBSULF_resampleconstant('ACTcarb',inputs,resampleN);

% rate ratio of chemical weathering in volcanic to non-volcanic silicate 
% rocks (called "Wv/Wnv in Berner 2006b, and "basalt/granite" in Berner 2008)
VNV = GEOCARBSULF_resampleconstant('VNV',inputs,resampleN);

% coefficient relating physical erosion to the mean 87Sr/86Sr of non-
% volcanic silicate rocks
NV = GEOCARBSULF_resampleconstant('NV',inputs,resampleN);

% exponent relating physical erosion to the mean 87Sr/86Sr of non-volcanic 
% silicate rocks (see Berner 2008)
exp_NV = GEOCARBSULF_resampleconstant('exp_NV',inputs,resampleN);

% rate ratio of chemical weathering in a minimally-vegetated to present-day 
% (angiosperm dominated) world
LIFE = GEOCARBSULF_resampleconstant('LIFE',inputs,resampleN);

% rate ratio of chemical weathering by gymnosperms to angiosperms
GYM = GEOCARBSULF_resampleconstant('GYM',inputs,resampleN);

% exponent reflecting the fraction of vegetation whose growth is stimulated 
% by elevated CO2; FERT is related to enhanced chemical weathering by the 
% Michaelis-Menton expression [2RCO2/(1+RCO2)]^FERT, which is called 
% fBb(CO2) in Berner 2004 (p. 24)
FERT = GEOCARBSULF_resampleconstant('FERT',inputs,resampleN);

% exponent used to describe the effect of climate on silicate or carbonate 
% weathering in the absence of vascular plants at time (t) relative to the 
% present-day (see pp. 25 & 53-54 in Berner 2004 and pp. 67-68 in Berner 1994)
exp_fnBb = GEOCARBSULF_resampleconstant('exp_fnBb',inputs,resampleN);

% climate sensitivity (K per CO2 doubling)
deltaT2X = GEOCARBSULF_resampleconstant('deltaT2X',inputs,resampleN);

% factor by which deltaT2X changes during times with large continental ice 
% sheets
GLAC = GEOCARBSULF_resampleconstant('GLAC',inputs,resampleN);

% coefficient used to calculate CAPd13C (called "alphac" in BASIC code), 
% the stable carbon isotopic fractionation between shallow-marine carbonate 
% and shallow-marine organic matter
J = GEOCARBSULF_resampleconstant('J',inputs,resampleN);

% coefficient used to calculate CAPd34S (called "alphas" in BASIC code), 
% the stable sulfur isotopic fractionation between marine sulfate sulfur 
% and marine pyrite sulfur
n = GEOCARBSULF_resampleconstant('n',inputs,resampleN);

% effect on temperature from the linear increase in solar luminosity over 
% time (K per 570 Myrs)
Ws = GEOCARBSULF_resampleconstant('Ws',inputs,resampleN);

% exponent that scales the dilution of dissolved HCO3- with runoff (fD) 
% (see pp. 29-31 & 34-36 in Berner 2004)
exp_fD = GEOCARBSULF_resampleconstant('exp_fD',inputs,resampleN);

%                   NOTE
% GEOCARBSULF is run forward in time; the following are 
% present-day (t = 0 Myrs ago) or initial values (t = 570 Myrs ago) of 
% parameters that are subsequently recalculated at each time-step 
% (see Berner 2004, 2006a for details)
  
%                   present-day values

% sulfate flux from oxidative weathering of old pyrite at present-day
Fwpa_0 = GEOCARBSULF_resampleconstant('Fwpa_0',inputs,resampleN);

% sulfate flux from weathering of CaSO4 sulfur at present-day
Fwsa_0 = GEOCARBSULF_resampleconstant('Fwsa_0',inputs,resampleN);

% carbon flux from weathering of old sedimentary organic matter at present-day
Fwga_0 = GEOCARBSULF_resampleconstant('Fwga_0',inputs,resampleN);

% carbon flux from weathering of old Ca and Mg carbonates at present-day
Fwca_0 = GEOCARBSULF_resampleconstant('Fwca_0',inputs,resampleN);

% carbon degassing flux from volcanism, metamorphism, and diagenesis of 
% organic matter at present-day
Fmg_0 = GEOCARBSULF_resampleconstant('Fmg_0',inputs,resampleN);

% carbon degassing flux from volcanism, metamorphism, and diagenesis of 
% carbonates at present-day
Fmc_0 = GEOCARBSULF_resampleconstant('Fmc_0',inputs,resampleN);

% sulfur degassing flux from volcanism, metamorphism, and diagenesis of 
% pyrite at present-day
Fmp_0 = GEOCARBSULF_resampleconstant('Fmp_0',inputs,resampleN);

% sulfur degassing flux from volcanism, metamorphism, and diagenesis of 
% CaSO4 sulfur at present-day
Fms_0 = GEOCARBSULF_resampleconstant('Fms_0',inputs,resampleN);

% weathering flux for all Ca and Mg silicates at present-day
Fwsi_0 = GEOCARBSULF_resampleconstant('Fwsi_0',inputs,resampleN);

% fraction of total Ca and Mg silicate weathering drived from volcanic 
% rocks at present-day
Xvolc_0 = GEOCARBSULF_resampleconstant('Xvolc_0',inputs,resampleN);

% stable carbon isotopic fractionation between shallow-marine carbonate and 
% shallow-marine organic matter at present-day
CAPd13C_0 = GEOCARBSULF_resampleconstant('CAPd13C_0',inputs,resampleN);

% stable sulfur isotopic fractionation between shallow-marine CaSO4 sulfur 
% and pyrite sulfur at present-day
CAPd34S_0 = GEOCARBSULF_resampleconstant('CAPd34S_0',inputs,resampleN);

%                   values at start time


% mass of atmospheric O2 at 570 Myrs ago (ROYER) or specified start time
oxygen_start = GEOCARBSULF_resampleconstant('oxygen_570',inputs,resampleN);

% mass of young crustal organic carbon at 570 Myrs ago
Gy_start = GEOCARBSULF_resampleconstant('Gy_570',inputs,resampleN);

% mass of young crustal carbonate carbon at 570 Myrs ago
Cy_start = GEOCARBSULF_resampleconstant('Cy_570',inputs,resampleN);

% mass of old crustal carbonate carbon at 570 Myrs ago
Ca_start = GEOCARBSULF_resampleconstant('Ca_570',inputs,resampleN);

% mass of young CaSO4 sulfur at 570 Myrs ago
Ssy_start = GEOCARBSULF_resampleconstant('Ssy_570',inputs,resampleN);

% mass of young pyrite sulfur at 570 Myrs ago
Spy_start = GEOCARBSULF_resampleconstant('Spy_570',inputs,resampleN);

% d34S of young CaSO4 sulfur at 570 Myrs ago
dlsy_start = GEOCARBSULF_resampleconstant('dlsy_570',inputs,resampleN);

% d13C of young carbonate carbon at 570 Myrs ago
dlcy_start = GEOCARBSULF_resampleconstant('dlcy_570',inputs,resampleN);

% d34S of young pyrite sulfur at 570 Myrs ago
dlpy_start = GEOCARBSULF_resampleconstant('dlpy_570',inputs,resampleN);

% d34S of old pyrite sulfur at 570 Myrs ago
dlpa_start = GEOCARBSULF_resampleconstant('dlpa_570',inputs,resampleN);

% d13C of young organic matter at 570 Myrs ago
dlgy_start = GEOCARBSULF_resampleconstant('dlgy_570',inputs,resampleN);

% d13C of old organic matter at 570 Myrs ago
dlga_start = GEOCARBSULF_resampleconstant('dlga_570',inputs,resampleN);

% 87Sr/86Sr of young carbonates undergoing weathering at 570 Myrs ago
Rcy_start = GEOCARBSULF_resampleconstant('Rcy_570',inputs,resampleN);

% 87Sr/86Sr of old carbonates undergoing weathering at 570 Myrs ago
Rca_start = GEOCARBSULF_resampleconstant('Rca_570',inputs,resampleN);

% 87Sr/86Sr of sub-aerial and submarine volcanic rocks at 570 Myrs ago
Rv_start = GEOCARBSULF_resampleconstant('Rv_570',inputs,resampleN);

% 87Sr/86Sr of non-volcanic silicates (granites) at 570 Myrs ago
Rg_start = GEOCARBSULF_resampleconstant('Rg_570',inputs,resampleN);

% OTHER 

% Ca and Mg flux between basalt and seawater
Fob =  GEOCARBSULF_resampleconstant('Fob',inputs,resampleN); 

% mass of carbon in ocean
COC = GEOCARBSULF_resampleconstant('COC',inputs,resampleN); 

% mass of old crustal organic carbon
Ga = GEOCARBSULF_resampleconstant('Ga',inputs,resampleN); 

% mass of old CaSO4 sulfur
Ssa = GEOCARBSULF_resampleconstant('Ssa',inputs,resampleN);

% mass of old pyrite sulfur
Spa = GEOCARBSULF_resampleconstant('Spa',inputs,resampleN);

% mass of sulfur in oceans + "interacting rocks" (i.e., carbon in rocks 
% undergoing weathering, burial, etc.)
ST = GEOCARBSULF_resampleconstant('ST',inputs,resampleN);

% d34S of ST
dlst = GEOCARBSULF_resampleconstant('dlst',inputs,resampleN);

% mass of carbon in oceans + "interacting rocks" (i.e., carbon in rocks 
% undergoing weathering, burial, etc.)
CT = GEOCARBSULF_resampleconstant('CT',inputs,resampleN);

% d13C of CT
dlct = GEOCARBSULF_resampleconstant('dlct',inputs,resampleN);

% rate constant expressing mass dependence for young pyrite sulfur
kwpy = GEOCARBSULF_resampleconstant('kwpy',inputs,resampleN);

% rate constant expressing mass dependence for young CaSO4 sulfur
kwsy = GEOCARBSULF_resampleconstant('kwsy',inputs,resampleN);

% rate constant expressing mass dependence for young organic matter weathering
kwgy = GEOCARBSULF_resampleconstant('kwgy',inputs,resampleN);

% rate constant expressing mass dependence for young carbonate weathering
kwcy = GEOCARBSULF_resampleconstant('kwcy',inputs,resampleN);

disp('Resampling and Filling of Input Arrays and Constants Complete.')











