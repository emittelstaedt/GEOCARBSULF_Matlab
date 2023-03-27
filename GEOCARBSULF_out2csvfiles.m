%% GEOCARBSULF - Set up and Save Output from simulations to csv files
%   Filename: GEOCARBSULF_out2csvfiles.m 
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
%  This script outputs here include the mean and percentiles for O2 and CO2 
% for each age plus all of the resampled values for the different 
% parameters so that the values used can be tracked. 
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
%% Process and Parse data Prior to output


%calculate median and percentiles for O2 and CO2
CO2 = median(CO2_resamples,2,'omitnan');
O2  = median( O2_resamples,2,'omitnan');

% percent failed runs
failed_runs = 100.*sum(isnan(CO2_resamples),2)./resampleN;

% Calculate percentile values 
% - values specified at top of GEOCARBSULF_main.m
if (resampleN>1)
    percentiles_CO2 = quantile(CO2_resamples, percentile_values, 2);
    percentiles_O2  = quantile( O2_resamples, percentile_values, 2);
end



%merge results and export summary file to working directory
if (resampleN==1)
    GEOCARB_output = [age(:), failed_runs(:), CO2(:), O2(:)];
else
    if (loop_parameters==TRUE)
        %start_column = y*h-(y-1);
        %end_column = y*h;
        %GEOCARB_output[,start_column:end_column] = cbind(age, failed_runs, CO2, percentiles_CO2, O2, percentiles_O2)
        %column_header[,start_column:end_column] = c("age (Myrs ago)","failed runs (%)",paste("CO2 for",input[h,"parameter"],"(ppm)",sep=" "),percentile_values,paste("O2 for",input[h,"parameter"],"(%)",sep=" "),percentile_values)

        % update user on where we are on loops
        %disp(['finished with parameter # ',num2str(h),' of ',num2str(z)]);
        disp('NOT CURRENTLY IMPLEMENTED - PARAMETER LOOPS!!')
    else
        GEOCARB_output = [age(:), failed_runs(:), CO2(:), percentiles_CO2(:),...
            O2(:), percentiles_O2(:)];
    end
end %if...else statement


%merge results and export summary file to working directory
if (resampleN==1)
    GEOCARB_output <- cbind(age, failed_runs, CO2, O2)
else
    if (loop_parameters==TRUE)
        %start_column = y*h-(y-1);
        %end_column = y*h;
        %GEOCARB_output(:,start_column:end_column) = [age(:), failed_runs(:),...
        %               CO2(:), percentiles_CO2(:), O2(:), percentiles_O2(:)];

        %column_header[,start_column:end_column] <- c("age (Myrs ago)",...
        %    "failed runs (%)",paste("CO2 for",input[h,"parameter"],"(ppm)",sep=" "),...
        %    percentile_values,paste("O2 for",input[h,"parameter"],"(%)",sep=" "),...
        %    percentile_values)

        % update user on where we are on loops
        %disp(['finished with parameter # ',num2str(h),' of ',num2str(z)]);

    else
        % output summary file
        GEOCARB_output = [age(:), failed_runs(:), CO2(:), percentiles_CO2(:),...
                          O2(:), percentiles_O2(:)];

    end
end % if...else statement


% HERE IS WHERE THE INPUT PARAMETER LOOP WILL END (ONCE IMPLEMENTED - h)

% summary output save
if strcmp(loop_parameters,'TRUE')  
    colnames(GEOCARB_output) = column_header 
else
    !!write.csv(GEOCARB_output, "GEOCARB_output.csv")
end


% OUTPUT!!!


if (input_distribution==TRUE & loop_parameters==FALSE) {
  resampled_input_constants <- matrix(nrow=length(input[,"parameter"]), ncol=2)
  colnames(resampled_input_constants) <- c("resampled_mean", "resampled_two_sigma")
  rownames(resampled_input_constants) <- input[, "parameter"]
  resampled_input_constants <- cbind(input[, -c(7:11)], resampled_input_constants)
  resampled_input_constants["ACT", "resampled_mean"] <- mean(ACT, na.rm=TRUE); resampled_input_constants["ACT", "resampled_two_sigma"] <- 2*sd(ACT, na.rm=TRUE)
  resampled_input_constants["ACTcarb", "resampled_mean"] <- mean(ACTcarb, na.rm=TRUE); resampled_input_constants["ACTcarb", "resampled_two_sigma"] <- 2*sd(ACTcarb, na.rm=TRUE)
  resampled_input_constants["VNV", "resampled_mean"] <- mean(VNV, na.rm=TRUE); resampled_input_constants["VNV", "resampled_two_sigma"] <- 2*sd(VNV, na.rm=TRUE)
  resampled_input_constants["NV", "resampled_mean"] <- mean(NV, na.rm=TRUE); resampled_input_constants["NV", "resampled_two_sigma"] <- 2*sd(NV, na.rm=TRUE)
  resampled_input_constants["exp_NV", "resampled_mean"] <- mean(exp_NV, na.rm=TRUE); resampled_input_constants["exp_NV", "resampled_two_sigma"] <- 2*sd(exp_NV, na.rm=TRUE)
  resampled_input_constants["LIFE", "resampled_mean"] <- mean(LIFE, na.rm=TRUE); resampled_input_constants["LIFE", "resampled_two_sigma"] <- 2*sd(LIFE, na.rm=TRUE)
  resampled_input_constants["GYM", "resampled_mean"] <- mean(GYM, na.rm=TRUE); resampled_input_constants["GYM", "resampled_two_sigma"] <- 2*sd(GYM, na.rm=TRUE)
  resampled_input_constants["FERT", "resampled_mean"] <- mean(FERT, na.rm=TRUE); resampled_input_constants["FERT", "resampled_two_sigma"] <- 2*sd(FERT, na.rm=TRUE)
  resampled_input_constants["exp_fnBb", "resampled_mean"] <- mean(exp_fnBb, na.rm=TRUE); resampled_input_constants["exp_fnBb", "resampled_two_sigma"] <- 2*sd(exp_fnBb, na.rm=TRUE)
  log_deltaT2X <- log(deltaT2X); resampled_input_constants["deltaT2X", "resampled_mean"] <- exp(mean(log_deltaT2X, na.rm=TRUE)); resampled_input_constants["deltaT2X", "resampled_two_sigma"] <- exp(2*sd(log_deltaT2X, na.rm=TRUE))
  resampled_input_constants["GLAC", "resampled_mean"] <- mean(GLAC, na.rm=TRUE); resampled_input_constants["GLAC", "resampled_two_sigma"] <- 2*sd(GLAC, na.rm=TRUE)
  resampled_input_constants["J", "resampled_mean"] <- mean(J, na.rm=TRUE); resampled_input_constants["J", "resampled_two_sigma"] <- 2*sd(J, na.rm=TRUE)
  resampled_input_constants["n", "resampled_mean"] <- mean(n, na.rm=TRUE); resampled_input_constants["n", "resampled_two_sigma"] <- 2*sd(n, na.rm=TRUE)
  resampled_input_constants["Ws", "resampled_mean"] <- mean(Ws, na.rm=TRUE); resampled_input_constants["Ws", "resampled_two_sigma"] <- 2*sd(Ws, na.rm=TRUE)
  resampled_input_constants["exp_fD", "resampled_mean"] <- mean(exp_fD, na.rm=TRUE); resampled_input_constants["exp_fD", "resampled_two_sigma"] <- 2*sd(exp_fD, na.rm=TRUE)
  resampled_input_constants["Fwpa_0", "resampled_mean"] <- mean(Fwpa_0, na.rm=TRUE); resampled_input_constants["Fwpa_0", "resampled_two_sigma"] <- 2*sd(Fwpa_0, na.rm=TRUE)
  resampled_input_constants["Fwsa_0", "resampled_mean"] <- mean(Fwsa_0, na.rm=TRUE); resampled_input_constants["Fwsa_0", "resampled_two_sigma"] <- 2*sd(Fwsa_0, na.rm=TRUE)
  resampled_input_constants["Fwga_0", "resampled_mean"] <- mean(Fwga_0, na.rm=TRUE); resampled_input_constants["Fwga_0", "resampled_two_sigma"] <- 2*sd(Fwga_0, na.rm=TRUE)
  resampled_input_constants["Fwca_0", "resampled_mean"] <- mean(Fwca_0, na.rm=TRUE); resampled_input_constants["Fwca_0", "resampled_two_sigma"] <- 2*sd(Fwca_0, na.rm=TRUE)
  resampled_input_constants["Fmg_0", "resampled_mean"] <- mean(Fmg_0, na.rm=TRUE); resampled_input_constants["Fmg_0", "resampled_two_sigma"] <- 2*sd(Fmg_0, na.rm=TRUE)
  resampled_input_constants["Fmc_0", "resampled_mean"] <- mean(Fmc_0, na.rm=TRUE); resampled_input_constants["Fmc_0", "resampled_two_sigma"] <- 2*sd(Fmc_0, na.rm=TRUE)
  resampled_input_constants["Fmp_0", "resampled_mean"] <- mean(Fmp_0, na.rm=TRUE); resampled_input_constants["Fmp_0", "resampled_two_sigma"] <- 2*sd(Fmp_0, na.rm=TRUE)
  resampled_input_constants["Fms_0", "resampled_mean"] <- mean(Fms_0, na.rm=TRUE); resampled_input_constants["Fms_0", "resampled_two_sigma"] <- 2*sd(Fms_0, na.rm=TRUE)
  resampled_input_constants["Fwsi_0", "resampled_mean"] <- mean(Fwsi_0, na.rm=TRUE); resampled_input_constants["Fwsi_0", "resampled_two_sigma"] <- 2*sd(Fwsi_0, na.rm=TRUE)
  resampled_input_constants["Xvolc_0", "resampled_mean"] <- mean(Xvolc_0, na.rm=TRUE); resampled_input_constants["Xvolc_0", "resampled_two_sigma"] <- 2*sd(Xvolc_0, na.rm=TRUE)
  resampled_input_constants["CAPd13C_0", "resampled_mean"] <- mean(CAPd13C_0, na.rm=TRUE); resampled_input_constants["CAPd13C_0", "resampled_two_sigma"] <- 2*sd(CAPd13C_0, na.rm=TRUE)
  resampled_input_constants["CAPd34S_0", "resampled_mean"] <- mean(CAPd34S_0, na.rm=TRUE); resampled_input_constants["CAPd34S_0", "resampled_two_sigma"] <- 2*sd(CAPd34S_0, na.rm=TRUE)
  resampled_input_constants["oxygen_570", "resampled_mean"] <- mean(oxygen_570, na.rm=TRUE); resampled_input_constants["oxygen_570", "resampled_two_sigma"] <- 2*sd(oxygen_570, na.rm=TRUE)
  resampled_input_constants["Gy_570", "resampled_mean"] <- mean(Gy_570, na.rm=TRUE); resampled_input_constants["Gy_570", "resampled_two_sigma"] <- 2*sd(Gy_570, na.rm=TRUE)
  resampled_input_constants["Cy_570", "resampled_mean"] <- mean(Cy_570, na.rm=TRUE); resampled_input_constants["Cy_570", "resampled_two_sigma"] <- 2*sd(Cy_570, na.rm=TRUE)
  resampled_input_constants["Ca_570", "resampled_mean"] <- mean(Ca_570, na.rm=TRUE); resampled_input_constants["Ca_570", "resampled_two_sigma"] <- 2*sd(Ca_570, na.rm=TRUE)
  resampled_input_constants["Ssy_570", "resampled_mean"] <- mean(Ssy_570, na.rm=TRUE); resampled_input_constants["Ssy_570", "resampled_two_sigma"] <- 2*sd(Ssy_570, na.rm=TRUE)
  resampled_input_constants["Spy_570", "resampled_mean"] <- mean(Spy_570, na.rm=TRUE); resampled_input_constants["Spy_570", "resampled_two_sigma"] <- 2*sd(Spy_570, na.rm=TRUE)
  resampled_input_constants["dlsy_570", "resampled_mean"] <- mean(dlsy_570, na.rm=TRUE); resampled_input_constants["dlsy_570", "resampled_two_sigma"] <- 2*sd(dlsy_570, na.rm=TRUE)
  resampled_input_constants["dlcy_570", "resampled_mean"] <- mean(dlcy_570, na.rm=TRUE); resampled_input_constants["dlcy_570", "resampled_two_sigma"] <- 2*sd(dlcy_570, na.rm=TRUE)
  resampled_input_constants["dlpy_570", "resampled_mean"] <- mean(dlpy_570, na.rm=TRUE); resampled_input_constants["dlpy_570", "resampled_two_sigma"] <- 2*sd(dlpy_570, na.rm=TRUE)
  resampled_input_constants["dlpa_570", "resampled_mean"] <- mean(dlpa_570, na.rm=TRUE); resampled_input_constants["dlpa_570", "resampled_two_sigma"] <- 2*sd(dlpa_570, na.rm=TRUE)
  resampled_input_constants["dlgy_570", "resampled_mean"] <- mean(dlgy_570, na.rm=TRUE); resampled_input_constants["dlgy_570", "resampled_two_sigma"] <- 2*sd(dlgy_570, na.rm=TRUE)
  resampled_input_constants["dlga_570", "resampled_mean"] <- mean(dlga_570, na.rm=TRUE); resampled_input_constants["dlga_570", "resampled_two_sigma"] <- 2*sd(dlga_570, na.rm=TRUE)
  resampled_input_constants["Rcy_570", "resampled_mean"] <- mean(Rcy_570, na.rm=TRUE); resampled_input_constants["Rcy_570", "resampled_two_sigma"] <- 2*sd(Rcy_570, na.rm=TRUE)
  resampled_input_constants["Rca_570", "resampled_mean"] <- mean(Rca_570, na.rm=TRUE); resampled_input_constants["Rca_570", "resampled_two_sigma"] <- 2*sd(Rca_570, na.rm=TRUE)
  resampled_input_constants["Rv_570", "resampled_mean"] <- mean(Rv_570, na.rm=TRUE); resampled_input_constants["Rv_570", "resampled_two_sigma"] <- 2*sd(Rv_570, na.rm=TRUE)
  resampled_input_constants["Rg_570", "resampled_mean"] <- mean(Rg_570, na.rm=TRUE); resampled_input_constants["Rg_570", "resampled_two_sigma"] <- 2*sd(Rg_570, na.rm=TRUE)
  resampled_input_constants["Fob", "resampled_mean"] <- mean(Fob, na.rm=TRUE); resampled_input_constants["Fob", "resampled_two_sigma"] <- 2*sd(Fob, na.rm=TRUE)
  resampled_input_constants["COC", "resampled_mean"] <- mean(COC, na.rm=TRUE); resampled_input_constants["COC", "resampled_two_sigma"] <- 2*sd(COC, na.rm=TRUE)
  resampled_input_constants["Ga", "resampled_mean"] <- mean(Ga, na.rm=TRUE); resampled_input_constants["Ga", "resampled_two_sigma"] <- 2*sd(Ga, na.rm=TRUE)
  resampled_input_constants["Ssa", "resampled_mean"] <- mean(Ssa, na.rm=TRUE); resampled_input_constants["Ssa", "resampled_two_sigma"] <- 2*sd(Ssa, na.rm=TRUE)
  resampled_input_constants["Spa", "resampled_mean"] <- mean(Spa, na.rm=TRUE); resampled_input_constants["Spa", "resampled_two_sigma"] <- 2*sd(Spa, na.rm=TRUE)
  resampled_input_constants["ST", "resampled_mean"] <- mean(ST, na.rm=TRUE); resampled_input_constants["ST", "resampled_two_sigma"] <- 2*sd(ST, na.rm=TRUE)
  resampled_input_constants["dlst", "resampled_mean"] <- mean(dlst, na.rm=TRUE); resampled_input_constants["dlst", "resampled_two_sigma"] <- 2*sd(dlst, na.rm=TRUE)
  resampled_input_constants["CT", "resampled_mean"] <- mean(CT, na.rm=TRUE); resampled_input_constants["CT", "resampled_two_sigma"] <- 2*sd(CT, na.rm=TRUE)
  resampled_input_constants["dlct", "resampled_mean"] <- mean(dlct, na.rm=TRUE); resampled_input_constants["dlct", "resampled_two_sigma"] <- 2*sd(dlct, na.rm=TRUE)
  resampled_input_constants["kwpy", "resampled_mean"] <- mean(kwpy, na.rm=TRUE); resampled_input_constants["kwpy", "resampled_two_sigma"] <- 2*sd(kwpy, na.rm=TRUE)
  resampled_input_constants["kwsy", "resampled_mean"] <- mean(kwsy, na.rm=TRUE); resampled_input_constants["kwsy", "resampled_two_sigma"] <- 2*sd(kwsy, na.rm=TRUE)
  resampled_input_constants["kwgy", "resampled_mean"] <- mean(kwgy, na.rm=TRUE); resampled_input_constants["kwgy", "resampled_two_sigma"] <- 2*sd(kwgy, na.rm=TRUE)
  resampled_input_constants["kwcy", "resampled_mean"] <- mean(kwcy, na.rm=TRUE); resampled_input_constants["kwcy", "resampled_two_sigma"] <- 2*sd(kwcy, na.rm=TRUE)
  
  resampled_input_arrays <- time_arrays
  for (j in 1:ageN) {
    resampled_input_arrays[j, match("Sr", colnames(resampled_input_arrays))] <- mean(Sr[j,], na.rm=TRUE); resampled_input_arrays[j, match("Sr", colnames(resampled_input_arrays))+1] <- 2*sd(Sr[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("d13C", colnames(resampled_input_arrays))] <- mean(d13C[j,], na.rm=TRUE); resampled_input_arrays[j, match("d13C", colnames(resampled_input_arrays))+1] <- 2*sd(d13C[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("d34S", colnames(resampled_input_arrays))] <- mean(d34S[j,], na.rm=TRUE); resampled_input_arrays[j, match("d34S", colnames(resampled_input_arrays))+1] <- 2*sd(d34S[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("fR", colnames(resampled_input_arrays))] <- mean(fR[j,], na.rm=TRUE); resampled_input_arrays[j, match("fR", colnames(resampled_input_arrays))+1] <- 2*sd(fR[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("fL", colnames(resampled_input_arrays))] <- mean(fL[j,], na.rm=TRUE); resampled_input_arrays[j, match("fL", colnames(resampled_input_arrays))+1] <- 2*sd(fL[j,], na.rm=TRUE)
    if (Godderis==TRUE) {
      resampled_input_arrays[j, match("fA_Godderis", colnames(resampled_input_arrays))] <- mean(fA[j,], na.rm=TRUE); resampled_input_arrays[j, match("fA_Godderis", colnames(resampled_input_arrays))+1] <- 2*sd(fA[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("fAw_fA_Godderis", colnames(resampled_input_arrays))] <- mean(fAw_fA[j,], na.rm=TRUE); resampled_input_arrays[j, match("fAw_fA_Godderis", colnames(resampled_input_arrays))+1] <- 2*sd(fAw_fA[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("fD_Godderis", colnames(resampled_input_arrays))] <- mean(fD[j,], na.rm=TRUE); resampled_input_arrays[j, match("fD_Godderis", colnames(resampled_input_arrays))+1] <- 2*sd(fD[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("GEOG_Godderis", colnames(resampled_input_arrays))] <- mean(GEOG[j,], na.rm=TRUE); resampled_input_arrays[j, match("GEOG_Godderis", colnames(resampled_input_arrays))+1] <- 2*sd(GEOG[j,], na.rm=TRUE) }
    else {
      resampled_input_arrays[j, match("fA", colnames(resampled_input_arrays))] <- mean(fA[j,], na.rm=TRUE); resampled_input_arrays[j, match("fA", colnames(resampled_input_arrays))+1] <- 2*sd(fA[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("fAw_fA", colnames(resampled_input_arrays))] <- mean(fAw_fA[j,], na.rm=TRUE); resampled_input_arrays[j, match("fAw_fA", colnames(resampled_input_arrays))+1] <- 2*sd(fAw_fA[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("fD", colnames(resampled_input_arrays))] <- mean(fD[j,], na.rm=TRUE); resampled_input_arrays[j, match("fD", colnames(resampled_input_arrays))+1] <- 2*sd(fD[j,], na.rm=TRUE)
      resampled_input_arrays[j, match("GEOG", colnames(resampled_input_arrays))] <- mean(GEOG[j,], na.rm=TRUE); resampled_input_arrays[j, match("GEOG", colnames(resampled_input_arrays))+1] <- 2*sd(GEOG[j,], na.rm=TRUE)
    }
    resampled_input_arrays[j, match("RT", colnames(resampled_input_arrays))] <- mean(RT[j,], na.rm=TRUE); resampled_input_arrays[j, match("RT", colnames(resampled_input_arrays))+1] <- 2*sd(RT[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("fSR", colnames(resampled_input_arrays))] <- mean(fSR[j,], na.rm=TRUE); resampled_input_arrays[j, match("fSR", colnames(resampled_input_arrays))+1] <- 2*sd(fSR[j,], na.rm=TRUE)
    resampled_input_arrays[j, match("fC", colnames(resampled_input_arrays))] <- mean(fC[j,], na.rm=TRUE); resampled_input_arrays[j, match("fC", colnames(resampled_input_arrays))+1] <- 2*sd(fC[j,], na.rm=TRUE)      
  }
  write.csv(resampled_input_arrays, "resampled_input_arrays.csv")
  write.csv(resampled_input_constants, "resampled_input_constants.csv")
}




