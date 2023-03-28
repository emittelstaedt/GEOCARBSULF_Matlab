%% GEOCARBSULF - Plot CO2 and O2 outputs from GEOCARBSULF run
%   Filename: GEOCARBSULF_summaryplots.m 
% 
%     See the GEOCARBSULF_main.m comments for notes on the calculations and
% implementation of GEOCARBSULF in matlab.  Here, the output from that code
% is taken and a summary plot of the atmospheric CO2 and O2 concentrations
% are plotted.  
%     If resampling was performed (resampleN>1), a confidence
% interval based upon the first two values of the "percentile_values" array  
% is plotted in the backgroud of the data.  For values of 0.025 and 0.975,
% this plot will be similar to Figure 2 in Royer et al., 2014. 
%
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
%% CO2 Plot

figure(1),clf;

subplot(2,1,1)

yconf = [percentiles_CO2(:,2); percentiles_CO2(end:-1:1,1)];
xconf = [age(:); age(end:-1:1)];
% plot confidence/uncertainty interval if resampleN>1
if (resampleN>1)
    p = fill(xconf,yconf,'black');
    p.FaceColor = [0.9 0.9 0.9];
    p.EdgeColor = 'none';
    hold on
end
semilogy(age,CO2,'k','LineWidth',2) 
set(gca,'Xdir','reverse')
set(gca,'Xlim',[0 570])
set(gca,'Yscale','log')
xlabel('Time (Ma)');
ylabel('Atmospheric CO_2 (ppm)')
title(['Results for ',num2str(irs-nfail),' successful runs.'])

subplot(2,1,2)

yconf = [percentiles_O2(:,2); percentiles_O2(end:-1:1,1)];
xconf = [age(:); age(end:-1:1)];
%set(gca,'')
% plot confidence/uncertainty interval if resampleN>1
if (resampleN>1)
    p = fill(xconf,yconf,'black');
    p.FaceColor = [0.9 0.9 0.9];
    p.EdgeColor = 'none';
    hold on
end
semilogy(age,O2,'k','LineWidth',2) 
set(gca,'Xdir','reverse')
set(gca,'Xlim',[0 570])
set(gca,'Yscale','log')
xlabel('Time (Ma)');
ylabel('Atmospheric O_2 (%)')





