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
plot_Royer = 'TRUE'; % plot Royer et al., 2014 outputs (autodetect resampleN>1 or not)
%% CO2 Plot

figure(1),clf;

subplot(2,1,1)

% plot confidence/uncertainty interval if resampleN>1
if (resampleN>1)
    yconf = [percentiles_CO2(:,2); percentiles_CO2(end:-1:1,1)];
    xconf = [age(:); age(end:-1:1)];

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
axis([0 570 100 1e5])

subplot(2,1,2)

% plot confidence/uncertainty interval if resampleN>1
if (resampleN>1)
    yconf = [percentiles_O2(:,2); percentiles_O2(end:-1:1,1)];
    xconf = [age(:); age(end:-1:1)];

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
axis([0 570 5 45])

%% ROYER OUTPUTS

if (strcmp(plot_Royer,'TRUE'))
    royer_out = 'D:\Dropbox\CODE_AND_SCRIPTS\GEOCARBSULF\GEOCARB_Royer2014\GEOCARB_output.csv';

    % Load Royer et al., 2014 R-code outputs
    routputs = readtable(royer_out);

    % resampleN > 1 in R code run
    if (width(routputs) > 5)
        % Parse resampleN>1 Royer Output
        Age_r2014 = routputs{2:end,2};
        CO2_r2014 = routputs{2:end,4};
        CO2_025_r2014 = routputs{2:end,5};
        CO2_975_r2014 = routputs{2:end,6};

        O2_r2014 = routputs{2:end,7};
        O2_025_r2014 = routputs{2:end,8};
        O2_975_r2014 = routputs{2:end,9};

        figure(2),clf;

        subplot(2,1,1)

        yconf = [CO2_025_r2014(:); CO2_975_r2014(end:-1:1)];
        xconf = [Age_r2014(:); Age_r2014(end:-1:1)];
        
        % plot confidence/uncertainty interval
        p = fill(xconf,yconf,'black');
        p.FaceColor = [0.9 0.9 0.9];
        p.EdgeColor = 'none';
        hold on

        semilogy(Age_r2014,CO2_r2014,'g','LineWidth',2)
        set(gca,'Xdir','reverse')
        set(gca,'Xlim',[0 570])
        set(gca,'Yscale','log')
        xlabel('Time (Ma)');
        ylabel('Atmospheric CO_2 (ppm)')
        title('Royer et al., 2014 R code - resampleN>1')
        axis([0 570 100 1e5])

        subplot(2,1,2)

        yconf = [O2_025_r2014(:); O2_975_r2014(end:-1:1)];
        xconf = [Age_r2014(:); Age_r2014(end:-1:1)];

        % plot confidence/uncertainty interval if resampleN>1
        p = fill(xconf,yconf,'black');
        p.FaceColor = [0.9 0.9 0.9];
        p.EdgeColor = 'none';
        hold on

        semilogy(Age_r2014,O2_r2014,'g','LineWidth',2)
        set(gca,'Xdir','reverse')
        set(gca,'Xlim',[0 570])
        set(gca,'Yscale','log')
        xlabel('Time (Ma)');
        ylabel('Atmospheric O_2 (%)')
        axis([0 570 5 45])

    % resample == 1
    else

        % Parse resampleN>1 Royer Output
        Age_r2014 = routputs{2:end,2};
        CO2_r2014 = routputs{2:end,4};
        O2_r2014 = routputs{2:end,5};

        figure(2),clf;

        subplot(2,1,1)
        semilogy(Age_r2014,CO2_r2014,'g','LineWidth',2)
        set(gca,'Xdir','reverse')
        set(gca,'Xlim',[0 570])
        set(gca,'Yscale','log')
        xlabel('Time (Ma)');
        ylabel('Atmospheric CO_2 (ppm)')
        title('Royer et al., 2014 R code - resampleN=1')
        axis([0 570 100 1e5])

        subplot(2,1,2)
        semilogy(Age_r2014,O2_r2014,'g','LineWidth',2)
        set(gca,'Xdir','reverse')
        set(gca,'Xlim',[0 570])
        set(gca,'Yscale','log')
        xlabel('Time (Ma)');
        ylabel('Atmospheric O_2 (%)')
        axis([0 570 5 45])
    end
end
