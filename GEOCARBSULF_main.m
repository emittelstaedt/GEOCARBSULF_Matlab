%% GEOCARBSULF w/Monte Carlo Resample Option 
%   Based upon the R script written by Dana Royer for the publication Royer 
% et al (2014). Please see the ROYER_README.txt file in the same directory
% as this script for details on their work and the R code. 
%
%   This script is a translation of the R code to Matlab with a
% reorganization into functions to make the overall code structure more
% readable and easier to understand than the single file R code.
% Additionally, a time-varying degassing term Xvolc was added to allow
% calculation of the CO2 contribution of hotspot volcanic fluxes.
%
% *******************    INPUT FILE DESCRIPTIONS    ***********************
%   The input files used in the R code are unchanged and are described by 
% Royer in their R code as follows: 
%   There are two input files: the default values reproduce the light gray
% uncertainty envelopes in Royer et al (2014); to reproduce the dark gray
% envelopes, reduce all two sigma errors by 75%.
%       "GEOCARB_input_summaries.csv" provides a summary of all input
%                                     parameters (see also Appendices 1-3 
%                                     in Royer et al, 2014), including
%                                     all necessary information for the 
%                                     time-invariant constants.
%       "GEOCARB_input_arrays.csv" provides the time-step mean and two sigma
%                                     values for the subset of parameters 
%                                     that are time arrays. 
% Please note:
%  If you open GEOCARB_input_summaries.csv in Excel, the "-inf" cells will
% probably be read in as an error; to correct, type <<'-inf>>.
% 
%   In the "input summaries" file: parameters ending in "_Godderis" are only
% activated if "Godderis" is set to TRUE (third line in the code below); in
% this case, the corresponding non-Godderis parameter choices are ignored.
%    The "type" column differentiates time-invariant parameters (constants)
%       from time-dependent parameters (arrays). 
%
%    The "resample" column describes whether a particular parameter will 
%       undergo resampling for error analysis. If you wish to forgo all 
%       resampling, set resampleN to 1 (first line in the code below); this 
%       will override whatever is written in this particular column. 
%
%    The "distribution" column describes whether the resampled distributions 
%       will follow a normal or lognormal distribution.
%
%    The "lower_limit" and "upper_limit" columns describe whether resamples
%       should be clipped at particular threshold values; in these cases, a 
%       value just inside the limit (+/- 0.0001) is assigned; in the default 
%       model run, this rarely happens.
%
% ******************    END INPUT FILE DESCRIPTIONS   ********************
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
%   For additional information including references to pertinent papers,
% please see the ROYER_README.txt file in this directory.
%
clear
%%  PARAMETERS

% ******************   MONTE CARLO ANALYSIS PARAMETERS   ******************
% number of resamples; if you wish to bypass resampling, set to 1
resampleN = 100;  

% list of percentiles (of any length) used to evaluate a resampled data set
% the median (0.5) is outputted by default 
percentile_values = [0.025,0.975];  

% set to "TRUE" to return the means and standard deviations of the input 
% parameter choices that are associated with successful (non-failed) runs; 
% only works when loop_parameters is set to "FALSE"; 
% this section of code was not used in Royer et al (2014) 
% (i.e., parameter set to "FALSE")
input_distribution = 'FALSE'; 

% set to "TRUE" to test the effect on calculated CO2 and O2 by sequentially 
% varying one input parameter at a time (figures 3-4 in Royer et al 2014); 
% BEWARE: long run times (68X longer than a single resampled run)
loop_parameters = 'FALSE'; 

% maximum number of times the convergence equation for CO2 will iterate 
% before signaling a failed run; in a test with reasonably-well-constrained 
% input parameters (similar to simulations presented in Royer et al, 2014), 
% the number of iterations never exceeded 7; this variable is here mostly 
% as a failsafe stop-gap.
iteration_threshold = 10; 

% ******************   GEOCARBSULF PROCESS PARAMETERS   *******************
% set to "TRUE" to run time arrays of fA, fAw/fA, fD, and GEOG from 
% Godderis et al, 2012; "FALSE" to run standard GEOCARBSULF time arrays 
Godderis = 'TRUE';  

% Mass of Present day atmospheric O2 (units???)
oxygen_0 = 38;

% Time step size (millions of years, Myrs) 
Dt = 10; 


%% Read Inputs and Prepare Output Files

% Load Input Files
GEOCARBSULF_loadinputfiles;


%% Setup Arrays for Calculations

% *********** Resample Arrays and Constants for MC error Analysis *********
% Note: if resampleN ==1 no resampling will be done
% if resample set to 'FALSE' in input files that particular parameter will
% not be resampled
GEOCARBSULF_performresampling;

% Initialize Empty calculation arrays for individual calcs of O2 and CO2
O2_resamples = zeros(nsteps, resampleN);
CO2_resamples = zeros(nsteps,resampleN);


%% GEOCARBSULF CALCULATIONS (mass balance and chemical reactions)

% Two loops here: 1) over resamples (i.e., number of model runs) and 2)
% over tsteps for each run.  Will see if I can vectorize one or the other
% (or both??).  

% loop over number of resamples
for irs = 1:resampleN

    % initialize simulation starting values (ROYER - 570 Myr)
    GEOCARBSULF_setinitialcondition;

    % loop over time steps
    for tt = 1:nsteps

        % set flag for failed run to false
        faild_run = 'FALSE';

        % get current time
        t = time_arrays(tt,"age").age;

        % set factors influenced by glacial vs. non-glacial state ("GCM")
        % glacial periods 260 - 330 Myrs ago and 35 - 0 Myrs ago;
        % BASIC code calls for a 270-340 Myrs ago interval, but see
        % Fielding et al 2008 (GSA Special Paper 441: 343-354) for justification
        if ((t<=330 && t>=260) || t< 36)
            % in GEOCARBSULF, deltaT2X = GCM*ln(2); capital gamma in Berner (2004)
            GCM = GLAC(irs).*deltaT2X(irs)/log(2);
        else
            GCM = deltaT2X(irs)/log(2);
        end

        % **********  calculate factors related to vegetation type  *******

        % vegetation = domination by non-vascular plants
        if (t<=570 && t>380)
            % effect of plants on weathering rate at time (t) to the present-day
            fE = LIFE(irs);

            %effect of CO2 on plant-assisted weathering for carbonates at
            % time (t) to the present-day
            fBB = (1+ACTcarb(irs).*GCM.*log(RCO2)-ACTcarb(irs).*Ws(irs).*(t/570)+...
                ACTcarb(irs).*GEOG(tt,irs)).*RCO2.^exp_fnBb(irs);

        elseif (t<=380 && t>350)
            %vegetation = ramp-up to gymnosperm domination
            fE = (GYM(irs)-LIFE(irs))*((380-t)/30)+LIFE(irs);
            fBB = ((380-t)/30).*((1 + ACTcarb(irs).*GCM.*log(RCO2)- ...
                ACTcarb(irs).*Ws(irs).*(t/570)+...
                ACTcarb(irs).*GEOG(tt,irs)).*(2*RCO2./(1+RCO2)).^FERT(irs))+...
                ((t-350)/30).*(1 + ACTcarb(irs).*GCM.*log(RCO2)-...
                ACTcarb(irs).*Ws(irs).*(t/570)+...
                ACTcarb(irs).*GEOG(tt,irs)).*RCO2.^exp_fnBb(irs);
        elseif (t<=350)
            fBB = (1+ACTcarb(irs)*GCM*log(RCO2)-ACTcarb(irs)*Ws(irs)*(t/570)+...
                ACTcarb(irs)*GEOG(tt,irs))*(2*RCO2/(1+RCO2))^FERT(irs);
        elseif (t<=350 && t>130)
            %vegetation = gymnosperm domination
            fE = GYM(irs);
        elseif (t<=130 && t>80)
            %vegetation = ramp-up to angiosperm domination
            fE = (1-GYM(irs))*((130-t)/50)+GYM(irs);
        elseif (t<=80)
            % vegetation = angiosperm domination
            fE = 1;
        end

        % ***************  CALCULATE SOURCE FLUXES ************************
        %(weathering and degassing); see pp. 5660-5661 in Berner (2006a) for
        % a discussion of the "Fwp", "Fws", and "Fwg" fluxes

        %weathering flux of young pyrite
        Fwpy = fA(tt,irs)*fR(tt,irs)*kwpy(irs)*Spy;

        %weathering flux of young CaSO4 sulfur
        Fwsy = fA(tt,irs)*fD(tt,irs)*kwsy(irs)*Ssy;

        %weathering flux of young organic carbon
        Fwgy = fA(tt,irs)*fR(tt,irs)*kwgy(irs)*Gy;

        %weathering flux of young carbonate carbon
        Fwcy = fA(tt,irs)*fD(tt,irs)*fL(tt,irs)*fE*fBB*kwcy(irs)*Cy;

        %weathering flux of old pyrite
        Fwpa = fR(tt,irs)*Fwpa_0(irs);

        %weathering flux of old CaSO4 sulfur
        Fwsa = fA(tt,irs)*fD(tt,irs)*Fwsa_0(irs);

        %weathering flux of old organic carbon
        Fwga = fR(tt,irs)*Fwga_0(irs);

        %weathering flux of old carbonate carbon
        Fwca = fA(tt,irs)*fD(tt,irs)*fL(tt,irs)*fE*fBB*Fwca_0(irs);

        %sulfur degassing flux from volcanism, metamorphism, and diagenesis of pyrite
        Fmp = fSR(tt,irs)*Fmp_0(irs);

        %sulfur degassing flux from volcanism, metamorphism, and diagenesis of CaSO4
        Fms = fSR(tt,irs)*Fms_0(irs);

        %carbon degassing flux from volcanism, metamorphism, and diagenesis of organics
        Fmg = fSR(tt,irs)*Fmg_0(irs);

        %carbon degassing flux from volcanism, metamorphism, and diagenesis of carbonate
        Fmc = fSR(tt,irs)*fC(tt,irs)*Fmc_0(irs);

        Fyop = Fwpa+Fmp;  %degassing + weathering flux of pyrite
        Fyos = Fwsa+Fms;  %degassing + weathering flux of CaSO4 sulfur
        Fyog = Fwga+Fmg;  %degassing + weathering flux of organic carbon
        Fyoc = Fwca+Fmc;  %degassing + weathering flux of carbonate carbon

        % *****************  CALCULATE SINK FLUXES ************************









    end

end








