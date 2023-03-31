%% GEOCARBSULF - perform iterative solution for CO2 concentration in atm 
%   Filename: GEOCARBSULF_iterativeCO2solve.m
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
% fwsi_no_climate (calculated outside of this script), RCO2 can be 
% determined; see pp. 72-76 in Berner (2004) for details. 
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
%% Start iterative solve for CO2

% start with doubled CO2 to ensure iterative solve does not coverge on the
% first step.
RCO2_old = 2*RCO2;

iteration_count = 0;

% if you want to plot convergence for a limited number of cases
%it_track = zeros(10,1);
%co2_track = it_track;

% Begin iteration with equations depending on the current time period (see
% comment in header above). Note that "570" here is not replaced with
% tstart because it is specific to the time period here, not to the time at
% which the user wants to start their simulation. 
if (t<=570 && t>380 && strcmp(failed_run,'FALSE'))
    while  (abs((RCO2/RCO2_old)-1) > 0.01)
        iteration_count = iteration_count+1;
        RCO2_old = RCO2;
        fwsi_climate = ( RCO2.^(exp_fnBb(irs) + ACT(irs)*GCM) ).*...
            ( (1 + RT(tt,irs).*GCM.*log(RCO2) -...
            RT(tt,irs).*Ws(irs)*(t/570) + ...
            RT(tt,irs).*GEOG(tt,irs) ).^exp_fD(irs) ).* ...
            exp(-ACT(irs).*Ws(irs).*(t/570)).*exp(ACT(irs).*GEOG(tt,irs));       

        W = ( (exp_fnBb(irs) + ACT(irs).*GCM).*...
            (RCO2.^(-exp_fnBb(irs) + ACT(irs)*GCM)) ).*...
            ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).*(t/570) + ...
            RT(tt,irs)*GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs).*Ws(irs).*...
            (t/570)).*exp(ACT(irs)*GEOG(tt,irs));

       V = ( RCO2.^(exp_fnBb(irs) + ACT(irs).*GCM) ).*exp_fD(irs).*( (1 + ...
            RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).*(t/570) + ...
            RT(tt,irs).*GEOG(tt,irs) ).^(-(1-exp_fD(irs))) ).*...
            (RT(tt,irs).*GCM./RCO2).*exp(-ACT(irs).*Ws(irs).*t/570).*...
            exp(ACT(irs)*GEOG(tt,irs));

        if (isnan(fwsi_climate+W+V) || iteration_count==iteration_threshold)
            failed_run = 'TRUE';
            break;
        end

        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V)))
            %damp the iteration to avoid overshoot
            RCO2 = RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V));  
        else
            %damp the iteration (convert the iteration to geometric 
            % shrinkage to avoid nonpositive value in overshoot)
            RCO2 = 0.2*RCO2;  
        end

    end %while loop
end %t>380 loop


if (t<=380 && t>350 && strcmp(failed_run,'FALSE'))  %the expressions for this time interval are more complex because the effects of plants on weathering are linearly mixed in; this helps to prevent model failure
    while  (abs(RCO2/RCO2_old-1) > 0.01)
        iteration_count = iteration_count+1;
        RCO2_old = RCO2;
        fwsi_climate_old = (RCO2.^(exp_fnBb(irs)+ACT(irs)*GCM)).*...
            ( (1 + RT(tt,irs).*GCM.*log(RCO2)-...
            RT(tt,irs).*Ws(irs)*(t/570) + ...
            RT(tt,irs).*GEOG(tt,irs)).^exp_fD(irs)).* ...
            exp(-ACT(irs).*Ws(irs).*(t/570)).*exp(ACT(irs).*GEOG(tt,irs));
                
        W_old = (exp_fnBb(irs) + ACT(irs).*GCM).*...
            (RCO2.^(-exp_fnBb(irs) + ACT(irs)*GCM) ).*...
            ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).*(t/570) + ...
            RT(tt,irs)*GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs).*Ws(irs).*...
            (t/570)).*exp(ACT(irs)*GEOG(tt,irs));

        V_old = (RCO2.^(exp_fnBb(irs) + ACT(irs).*GCM)).*exp_fD(irs).*((1 + ...
            RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).*(t/570) + ...
            RT(tt,irs).*GEOG(tt,irs)).^(-(1-exp_fD(irs))) ).*...
            (RT(tt,irs).*GCM./RCO2).*exp(-ACT(irs).*Ws(irs).*(t/570)).*...
            exp(ACT(irs)*GEOG(tt,irs));


        fwsi_climate_new = ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ...
                           ACT(irs).*GCM)) ).*((1 + RCO2).^(-FERT(irs))).*...
                           ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).* ...
                           Ws(irs).*(t/570) + RT(tt,irs).* ...
                           GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs)* ...
                           Ws(irs).*(t/570)).*exp(ACT(irs)*GEOG(tt,irs));

        W_new = (2.^FERT(irs)).*(FERT(irs) + ACT(irs).*GCM)*...
                (RCO2.^(FERT(irs) + ACT(irs).*GCM - 1)).*...
                ((1 + RCO2).^-FERT(irs)).*((1 + RT(tt,irs).*GCM.*log(RCO2) - ...
                RT(tt,irs)*Ws(irs).*(t/570) + RT(tt,irs).* ...
                GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs).*Ws(irs).*...
                (t/570)).*exp(ACT(irs).*GEOG(tt,irs));

        V_new = ( -FERT(irs).*((1 + RCO2).^(-(1 + FERT(irs)))) ).* ...
                ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ACT(irs).*GCM)) ).* ...
                ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).* ...
                (t/570) + RT(tt,irs).*GEOG(tt,irs)).^exp_fD(irs)).* ...
                exp(-ACT(irs).*Ws(irs).*(t/570)).*exp(ACT(irs).*GEOG(tt,irs));

        X_new = exp_fD(irs).*((1 + RT(tt,irs).*GCM.*log(RCO2) - ...
                RT(tt,irs).*Ws(irs).*(t/570) + RT(tt,irs).*GEOG(tt,irs)).^ ...
                (-(1 - exp_fD(irs))) ).*(RT(tt,irs).*GCM./RCO2).* ...
                ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ACT(irs).*GCM)) ).* ...
                ((1 + RCO2).^-FERT(irs)).*exp(-ACT(irs).*Ws(irs).*(t/570)).* ...
                exp(ACT(irs).*GEOG(tt,irs));

        fwsi_climate = ((t - 350)./31).*fwsi_climate_old + ...
                       ((381 - t)./31).*fwsi_climate_new;

        Fw_v_x = ((t - 350)./31).*(W_old + V_old) + ...
                 ((381 - t)./31).*(W_new + V_new + X_new);


        % temporary failed test variable
        tmp = fwsi_climate_old + W_old + V_old + fwsi_climate_new + W_new +...
              V_new + X_new + fwsi_climate + Fw_v_x;

        if (isnan(tmp) || iteration_count==iteration_threshold)  
            failed_run = 'TRUE';
            break;
        end

        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V)))
            %damp the iteration to avoid overshoot
            RCO2 = RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V));  
        else
            %damp the iteration (convert the iteration to geometric 
            % shrinkage to avoid nonpositive value in overshoot)
            RCO2 = 0.2*RCO2;  
        end

    end %while loop
end % t<=380 & t>350 loop


if (t<=350 && strcmp(failed_run,'FALSE')) 
    while  (abs(RCO2/RCO2_old-1) > 0.01) 
        iteration_count = iteration_count+1;
        RCO2_old = RCO2;

        fwsi_climate = ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ...
            ACT(irs).*GCM)) ).*((1 + RCO2).^(-FERT(irs))).*...
            ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).* ...
            Ws(irs).*(t/570) + RT(tt,irs).* ...
            GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs)* ...
            Ws(irs).*(t/570)).*exp(ACT(irs)*GEOG(tt,irs));        
 
        W = (2.^FERT(irs)).*(FERT(irs) + ACT(irs).*GCM)*...
                (RCO2.^(FERT(irs) + ACT(irs).*GCM - 1)).*...
                ((1 + RCO2).^-FERT(irs)).*((1 + RT(tt,irs).*GCM.*log(RCO2) - ...
                RT(tt,irs)*Ws(irs).*(t/570) + RT(tt,irs).* ...
                GEOG(tt,irs)).^exp_fD(irs)).*exp(-ACT(irs).*Ws(irs).*...
                (t/570)).*exp(ACT(irs).*GEOG(tt,irs));
       
        V = ( -FERT(irs).*((1 + RCO2).^-(1 + FERT(irs))) ).* ...
                ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ACT(irs).*GCM)) ).* ...
                ((1 + RT(tt,irs).*GCM.*log(RCO2) - RT(tt,irs).*Ws(irs).* ...
                (t/570) + RT(tt,irs).*GEOG(tt,irs)).^exp_fD(irs)).* ...
                exp(-ACT(irs).*Ws(irs).*(t/570)).*exp(ACT(irs).*GEOG(tt,irs));

        X = exp_fD(irs).*( (1 + RT(tt,irs).*GCM.*log(RCO2) - ...
                RT(tt,irs).*Ws(irs).*(t/570) + RT(tt,irs).*GEOG(tt,irs)).^ ...
                -(1 - exp_fD(irs)) ).*(RT(tt,irs).*GCM./RCO2).* ...
                ( (2.^FERT(irs)).*(RCO2.^(FERT(irs) + ACT(irs).*GCM)) ).* ...
                ((1 + RCO2).^-FERT(irs)).*exp(-ACT(irs).*Ws(irs).*(t/570)).* ...
                exp(ACT(irs).*GEOG(tt,irs));        
       
        % Check for failed run
        if (isnan(fwsi_climate + W + V + X) || iteration_count==iteration_threshold)
            failed_run = 'TRUE';
            %test_plot_convergence;
            break;
        end

        if (RCO2 > ((fwsi_climate - fwsi_no_climate)/(W+V)))
            %damp the iteration to avoid overshoot
            RCO2 = RCO2-0.9*((fwsi_climate - fwsi_no_climate)/(W+V));  
        else
            %damp the iteration (convert the iteration to geometric 
            % shrinkage to avoid nonpositive value in overshoot)
            RCO2 = 0.2*RCO2;  
        end

        % For plotting convergence 
        %RCO2_ratio = RCO2./RCO2_old - 1;
        %it_track(iteration_count) = RCO2_ratio;
        %co2_track(iteration_count) = RCO2;
        %if (abs((RCO2/RCO2_old)-1) <= 0.01)
        %    test_plot_convergence;
        %end


    end % while loop
end % t<=350 loop

% END OF CO2 ITERATION SOLUTION


%test for failed runs => CO2 < 150 ppm or > 50000 ppm
if (RCO2<0.6 || RCO2>200)  
    failed_run = 'TRUE';
end









