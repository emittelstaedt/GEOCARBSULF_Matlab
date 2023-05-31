%% Change Timesteps on Inputs to GEOCARBSULF from Royer et al., 2014 
%   GEOCARBSULF_InputTimestepMOD.m
%   
%   This script reads in the input parameters in the csv files provided by
%   Royer et al., 2014 and interpolates to a finer timestep size.  This
%   will mean that the assumption of steady state for each timestep will be
%   worse for sulphur, but we are interested in CO2.  As long as we are
%   looking at several million years, the hope is that steady state will
%   still be a valid assumption. In process is a literature search to
%   justify a shorter timestep for CO2 alone. 
%
%   Also an option of this script is to add in the hotspot magma degassing
%   flux as calculated for Mittelstaedt et al., 2023.  THis covers only the
%   last 80 Myr (limited by hotspot chain lengths/ages) and much of the
%   data comes from Morrow and Mittelstaedt, 2021.  However, several other
%   estimates are also added for Iceland, the Columbia River Basalts, and
%   the Ethiopian (Afar) flood basalts (see CarbonInPlumes_AddNonMM21Data.m 
%   for references and data sources. If this data is desired, it will be 
%   decimated to the desired timestep for GEOCARBSULF (<250kyr not 
%   not recommended since hotspot data sampled at 50kyr). 
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
%   In the "input summaries" file: 
%    The "type" column differentiates time-invariant parameters (constants)
%       from time-dependent parameters (arrays). 
%
%    The "resample" column describes whether a particular parameter will 
%       undergo resampling for error analysis. 
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
%   For additional information including references to pertinent papers,
% please see the ROYER_README.txt file in this directory.
%
clear
%%  Parameters 

% Choose Desired Time Range and Timestep (units: Millions of Years)
tstart_desired = 570;
Dt_desired = 1; 


%% Load Input Files

% set to "FALSE" - otherwise load script will modify values (see
% GEOCARBSULF_main.m for details. DO NOT CHANGE HERE!
loop_parameters = 'FALSE'; 

% Load Input Files
GEOCARBSULF_loadinputfiles;

% Hotspot data switches 
add_hs_only = 0; %no interp, only hs data add
add_HSdegas = 'yes';  % add in hotspot degassing as an array 'yes', 'no'
hsdegas_table = 'GEOCARB_input_hsdegas.csv';

%% Grab current values from existing time-arrays

% set time parameters from inputs
% Time step size (millions of years, Myrs) - ASSUMED CONSTANT!!!! 
Dt = time_arrays(1,"age").age - time_arrays(2,"age").age; 

% start time (Myrs)
tstart = time_arrays(1,"age").age;

% check if need to do anything 
if ( (Dt == Dt_desired) && (tstart == tstart_desired))    
    if (width(time_arrays) == 33) 
        if strcmp(add_HSdegas,'yes')
            disp('Hotspot degassing fluxes not added.  Adding fluxes.')
            add_hs_only = 1;
        else
            error('Input files match desired times. No hotspot data.');
        end
    else
       error('Input files match desired times and already contain hotspot degassing data.');
    end
end

%% Arrange data and Interpolate

% extract existing time vector
told = time_arrays(:,"age").age;

% create new time vector
tnew = [tstart_desired:-Dt_desired:0];

% extract data from table format to an array
GCSinputs = table2array(time_arrays);
GCSnames = time_arrays.Properties.VariableNames;

% establish outputs array
OutputArray = zeros(length(tnew),width(time_arrays));

% input new time values
OutputArray(:,1) = tnew(:);

% interpolate to new times 
for jj = 2:width(time_arrays)
     OutputArray(:,jj) = interp1(told(:),GCSinputs(:,jj),tnew,'linear');

end

%% Add in Hotspot Degassing Data 

if strcmp(add_HSdegas,'yes')
    dt_gcs = Dt_desired;
    inputs = readtable(hsdegas_table);

    % extract data from table format to an array
    HSinputs = table2array(inputs);
    HSnames = inputs.Properties.VariableNames;

    % set up array space (volc only past 80Ma, but GEOCARBSULF data starts
    % at 570 Ma
    OutputArray(:,jj+1:jj+4) = zeros(length(tnew),4);

    % put hotspot degassing data into the output array for saving to
    % GEOCARBSULF csv file
    OutputArray(end-(height(inputs)-1):end,jj+1:jj+4) = HSinputs;

    % concat names to GCS names cell array
    GCSnames = [GCSnames,HSnames];

end


%% Save Data into New Table

GCSoutTable = array2table(OutputArray,'VariableNames',GCSnames);

array_out = cast(array_table,'char');
fname = [array_out(1:end-4),'_tMod.csv'];
writetable(GCSoutTable,fname);


%% TESTING

test_var = 0; 

if (test_var > 1)
    figure(10),clf;
    plot(told,GCSinputs(:,test_var),'r')
    hold on
    plot(told,GCSinputs(:,test_var),'rs')
    plot(OutputArray(:,1),OutputArray(:,test_var),'g')
    plot(OutputArray(:,1),OutputArray(:,test_var),'g.')
    xlabel('Age (Myr)');
    ylabel(GCSnames{test_var});

end





