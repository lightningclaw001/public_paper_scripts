clear;clc;

maxNumCompThreads(1);

% invert HPPC data

% ??????????????? data loaded format: c data at each discretization, for
% 1:200; chemical potential; reaction rate; indicator pulse data: +i, -i is
% indicator relaxation data: +i+N, -i-N is indicator

[xPulse, xRelax] = invertHPPC();




function [paramsOptim, dataRef, params_nondim, params] = loadDataMPETSimsFullCell()
%%% load data from single particle simulations
% returns: loaded current and voltage data for current pulses and
% relaxation segments

%R0_bar = 0.0007350654435686459 ish;


% pull one of the datasets...any of them
load('HPPC_data_5meV_pulse.mat');

% set voltage params of interest
voltRange = [0.1];
%voltRange = [0.01, 0.05, 0.1, 0.5];

%% Define Key Dimensions
% P is number of parameters of interest.
P = 6;
% N is number of cycles for each simulation.
N = params.num_cycles;
% V is number of voltage pulses in the HPPC cycle.
xVar = size(voltRange, 2);
% L is number of discretizations in particle
L = params.m;

% define small parameter
eps0 = sqrt(eps);

% set voltage range for each parameter
voltage_pulse_range = cell(1, N);
voltage_pulse_range(:) = {voltRange};
x_charge = cell(1, N);
x_charge(:) = {0};

paramsOptim = struct('voltage_pulse_range', voltage_pulse_range);
dataRef = struct('x_charge', x_charge);

%reference solutions
ctilde = 1:-0.0111:0.9;
% ctilde = 1:-0.0222:0.8;
clyte = 1000:-11.1:900;
Rf = 0:0.025:0.26;
% Rf = 0:0.025:0.26;

load('HPPC_data_5meV_pulse.mat');

kToe = 0.0257;
params_nondim.kToe = kToe;

% input these values from MPET simulation
L_c = 50e-6; L_a = 87.8e-6; poros_c = 0.4; poros_a = 0.4;
P_L_c = 0.69; P_L_a = 0.69;

Crate_redim = 12703.292014982677;
mpet_t_ref = 7.792736808950305;

%let's save the t_ref and average particle size
params_nondim.mpet_t_ref = mpet_t_ref;
particle_sizes = [1.02042229e-07, 1.05002089e-07, 9.06185587e-08, 8.61786715e-08, ...
    9.49947598e-08, 9.27263684e-08, 8.68791247e-08, 1.07055739e-07, 9.09937621e-08, ...
    1.00786570e-07, 8.99585960e-08, 1.11841919e-07, 9.86287689e-08, 9.58665055e-08, ...
    9.52893031e-08, 1.05600255e-07, 1.08275568e-07, 1.12523738e-07, 1.08473145e-07, 1.12427833e-07];
params_nondim.avg_part_size = mean(particle_sizes)/1e-7;
%params_nondim.avg_part_size = average_particle_size(sizes);
params_nondim.electrode = "c";

params_nondim.rescale_c = L_c * (1-poros_c)*P_L_c*mean(3./particle_sizes);
%for full cell
params_nondim.rescale_a = L_a * (1-poros_a)*P_L_a*mean(3./particle_sizes);params_nondim.ca_ratio = params_nondim.rescale_c/params_nondim.rescale_a;
params_nondim.electrolyte = 'SM';
%params_nondim.electrolyte = 'dilute';

for i1 = 1:xVar 
    % voltage pulse size
       
    % guess value at each state degradation
    for i2 = 1:10
%degradation state
        % load data for all sols for i2-th state of degradation x is time y
        % is state variable

%          load('C:\Users\ChemeGrad2019\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\dezhuang\projects\hppc_data\MPET_full_cell_current_pulse\all_deg_sphere\current_pulses_'+ ...
%             string(paramsOptim(1).voltage_pulse_range(i1)) + ...
%             'V\ctilde_' + sprintf('%.4f',ctilde(i2))+'_clyte_' + ...
%             sprintf('%.1f',clyte(i2)) + '_Rf_'+ sprintf('%.3f',Rf(i2)) ...
%             + '\sim_output\output_data.mat');

         load('C:\Users\ChemeGrad2019\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\dezhuang\projects\hppc_data\MPET_full_cell_current_pulse\c_lyte_only_sphere\current_pulses_'+ ...
            string(paramsOptim(1).voltage_pulse_range(i1)) + ...
            'V\clyte_' + ...
            sprintf('%.1f',clyte(i2))  ...
            + '\sim_output\output_data.mat');

         % load('C:\Users\ChemeGrad2019\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\dezhuang\projects\hppc_data\MPET_full_cell_current_pulse\R_f_a_only_sphere\current_pulses_'+ ...
         %    string(paramsOptim(1).voltage_pulse_range(i1)) + ...
         %    'V\Rf_'+ sprintf('%.3f',Rf(i2)) + '\sim_output\output_data.mat');
         % 

         % load('C:\Users\ChemeGrad2019\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\dezhuang\projects\hppc_data\MPET_full_cell_current_pulse\c_tilde_c_only_sphere\current_pulses_'+ ...
         %    string(paramsOptim(1).voltage_pulse_range(i1)) + ...
         %    'V\ctilde_' + sprintf('%.4f',ctilde(i2))+'\sim_output\output_data.mat');


        if i1 == 1 & i2 == 1
            % only get the reference values at the first state
            params_nondim.electrode = "c";
            params_nondim.muR_ref = 0; params_nondim.muR_ref_a = 0; params_nondim.muR_ref_c = 0;
            params_nondim.muR_ref_c = -mu_c(ffrac_c(1), params_nondim);
            params_nondim.electrode = "a"; params_nondim.muR_ref_a = -mu_c(ffrac_a(1), params_nondim);
            params_nondim.muR_ref = params_nondim.muR_ref_c - params_nondim.muR_ref_a;
        end 

        %SM models: phi_lyte does not electrolyte concentration
        %dilute models: phi_lyte includes electrolyte activity

        % iterate over all cycles
        for i3 = 1:N    
            %state of charge

            % pull charge data last column of y array contains charge state
            % information. charging is at to +i3, charge relaxation is at
            % i3+N. we also remove the indexes where there is no current
            idxCharge = (CCCVCPcycle_maccor_cycle_counter == 2+(i3-1)*4 & (~nearZero(current, 1e-16)));
            % save charge pulse time and data
            xChargePulse = phi_applied_times(idxCharge)*mpet_t_ref;
            % let's also check for repeated values
            check_equality = logical([(xChargePulse(2:end)-xChargePulse(1:end-1) == 0), 0]);
            %save current data redeminsionalized
            % current in values of A/m^2
            Rxn_charge_pulse = current(idxCharge) * Crate_redim;
            %save voltage data redimensionalized
            mures_charge_pulse = phi_applied(idxCharge);
            %for dilute model
            
            mures_charge_pulse_c = -mean(phi_bulk_c(idxCharge,:),2) + mean(phi_lyte_c(idxCharge,:),2);
            mures_charge_pulse_a = -mean(phi_bulk_a(idxCharge,:),2) + mean(phi_lyte_a(idxCharge,:),2);
            %for concentrated model
            if params_nondim.electrolyte == "SM"
                mures_charge_pulse_c = -mean(phi_bulk_c(idxCharge,:),2) + mean(phi_lyte_c(idxCharge,:),2) - log(a_l(mean(c_lyte_c(idxCharge,:),2), params_nondim));
                mures_charge_pulse_a = -mean(phi_bulk_a(idxCharge,:),2) + mean(phi_lyte_a(idxCharge,:),2) - log(a_l(mean(c_lyte_a(idxCharge,:),2), params_nondim));
            end

            %save SOC data
            c_charge_c = ffrac_c(idxCharge);
            c_charge_a = ffrac_a(idxCharge);

            %remove repeated indices]
            xChargePulse(check_equality) = [];
            Rxn_charge_pulse(check_equality) = [];
            mures_charge_pulse(check_equality) = [];
            mures_charge_pulse_c(check_equality) = [];
            mures_charge_pulse_a(check_equality) = [];

            %save relaxation charge data
            idxCharge = (CCCVCPcycle_maccor_cycle_counter == 3+(i3-1)*4);
            mures_charge_rlx = phi_applied(idxCharge);
            % for dilute model
            mures_charge_rlx_c = -mean(phi_bulk_c(idxCharge,:),2) + mean(phi_lyte_c(idxCharge,:),2);
            mures_charge_rlx_a = -mean(phi_bulk_a(idxCharge,:),2) + mean(phi_lyte_a(idxCharge,:),2);
            % for concentrated model
            if params_nondim.electrolyte == "SM"
                mures_charge_rlx_c = -mean(phi_bulk_c(idxCharge,:),2) + mean(phi_lyte_c(idxCharge,:),2) - log(a_l(mean(c_lyte_c(idxCharge,:),2), params_nondim));
                mures_charge_rlx_a = -mean(phi_bulk_a(idxCharge,:),2) + mean(phi_lyte_a(idxCharge,:),2) - log(a_l(mean(c_lyte_a(idxCharge,:),2), params_nondim));
            end 
             
            % pull discharge data last column of y array contains charge
            % state information. discharging is at to -i3, discharge
            % relaxation is at -i3-N. save discharge pulse time and data we
            % also remove the indexes where there is no current
            idxDischarge = (CCCVCPcycle_maccor_cycle_counter == i3*4 & (~nearZero(current, 1e-16)));
            % save charge pulse time and data
            xDischargePulse = phi_applied_times(idxDischarge)*mpet_t_ref;
            % let's also check for repeated values
            check_equality = logical([(xDischargePulse(2:end)-xDischargePulse(1:end-1) == 0), 0]);
            %save current data redeminsionalized
            % current in values of A/m^2
            Rxn_discharge_pulse = current(idxDischarge) * Crate_redim;

            %save voltage data redimensionalized
            mures_discharge_pulse = phi_applied(idxDischarge);
            % for dilute model
            mures_discharge_pulse_c = -mean(phi_bulk_c(idxDischarge,:),2) + mean(phi_lyte_c(idxDischarge,:),2);
            mures_discharge_pulse_a = -mean(phi_bulk_a(idxDischarge,:),2) + mean(phi_lyte_a(idxDischarge,:),2);
            % for concentrated model
            if params_nondim.electrolyte == "SM"
                mures_discharge_pulse_c = -mean(phi_bulk_c(idxDischarge,:),2) + mean(phi_lyte_c(idxDischarge,:),2) - log(a_l(mean(c_lyte_c(idxDischarge,:),2), params_nondim));
                mures_discharge_pulse_a = -mean(phi_bulk_a(idxDischarge,:),2) + mean(phi_lyte_a(idxDischarge,:),2) - log(a_l(mean(c_lyte_a(idxDischarge,:),2), params_nondim));
            end 
            %save SOC value
            c_discharge_c = ffrac_c(idxDischarge);
            c_discharge_a = ffrac_a(idxDischarge);                

            %remove repeated indices
            xDischargePulse(check_equality) = [];
            Rxn_discharge_pulse(check_equality) = [];
            mures_discharge_pulse(check_equality) = [];
            mures_discharge_pulse_c(check_equality) = [];
            mures_discharge_pulse_a(check_equality) = [];

            %save relaxation charge data
            idxDischarge = (CCCVCPcycle_maccor_cycle_counter == 4*i3+1);
            mures_discharge_rlx = phi_applied(idxDischarge);
            % for dilute model
            mures_discharge_rlx_c = -mean(phi_bulk_c(idxDischarge,:),2) + mean(phi_lyte_c(idxDischarge,:),2);
            mures_discharge_rlx_a = -mean(phi_bulk_a(idxDischarge,:),2) + mean(phi_lyte_a(idxDischarge,:),2);
            % for concentrated model
            if params_nondim.electrolyte == "SM"
                mures_discharge_rlx_c = -mean(phi_bulk_c(idxDischarge,:),2) + mean(phi_lyte_c(idxDischarge,:),2) - log(a_l(mean(c_lyte_c(idxDischarge,:),2), params_nondim));
                mures_discharge_rlx_a = -mean(phi_bulk_a(idxDischarge,:),2) + mean(phi_lyte_a(idxDischarge,:),2) - log(a_l(mean(c_lyte_a(idxDischarge,:),2), params_nondim));
            end 

            if i1 == 1 && i2 == 1
                %preallocate arrays for
                dataRef(i3).y_charge_pulse = zeros(length(paramsOptim(1).voltage_pulse_range),...
                    10,length(xChargePulse));
                dataRef(i3).y_discharge_pulse = zeros(length(paramsOptim(1).voltage_pulse_range),...
                    10,length(xDischargePulse));

                %reference solution for discharge relaxation voltage
                %difference (only one data point at t =0)
                dataRef(i3).discharge_rlx = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                
                paramsOptim(i3).t_charge = xChargePulse-xChargePulse(1)+1e-8;
                paramsOptim(i3).t_discharge = (xDischargePulse-xDischargePulse(1) + 1e-8);

                %optimization parameters used--charge state for charge and
                %discharge curves
                paramsOptim(i3).c_charge_c = zeros(length(paramsOptim(1).voltage_pulse_range),2);
                paramsOptim(i3).c_discharge_c = zeros(length(paramsOptim(1).voltage_pulse_range),2);
                paramsOptim(i3).c_charge_a = zeros(length(paramsOptim(1).voltage_pulse_range),2);
                paramsOptim(i3).c_discharge_a = zeros(length(paramsOptim(1).voltage_pulse_range),2);
                paramsOptim(i3).charge_rlx = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).discharge_rlx = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_charge = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
                paramsOptim(i3).mures_discharge = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
                paramsOptim(i3).mures_charge_c = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
                paramsOptim(i3).mures_charge_a = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
                paramsOptim(i3).mures_discharge_c = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
                paramsOptim(i3).mures_discharge_a = zeros(length(paramsOptim(1).voltage_pulse_range),length(xChargePulse));
           end                 

            if i2 == 1
                %when the system has not degraded yet store the known
                %parameters
                paramsOptim(i3).c_charge_c(i1,1) = c_charge_c(1);
                paramsOptim(i3).c_discharge_c(i1,1) = c_discharge_c(1);
                paramsOptim(i3).c_charge_a(i1,1) = c_charge_a(1);
                paramsOptim(i3).c_discharge_a(i1,1) = c_discharge_a(1);                
                %after the charge/discharge pulse
                paramsOptim(i3).c_charge_c(i1,2) = c_charge_c(end);
                paramsOptim(i3).c_discharge_c(i1,2) = c_discharge_c(end);
                paramsOptim(i3).c_charge_a(i1,2) = c_charge_a(end);
                paramsOptim(i3).c_discharge_a(i1,2) = c_discharge_a(end);    
                %before the charge/discharge pulse
             
                paramsOptim(i3).c_charge_c(i1,1) = c_charge_c(1);
                paramsOptim(i3).c_discharge_c(i1,1) = c_discharge_c(1);          
                %after the charge/discharge pulse
                paramsOptim(i3).c_charge_c(i1,2) = c_charge_c(end);
                paramsOptim(i3).c_discharge_c(i1,2) = c_discharge_c(end);  
                %before the charge/discharge pulse
             
                %store the time data for charge/discharge in the
                %optim_params
               
                paramsOptim(i3).Rxn_charge(i1) = Rxn_charge_pulse(2); %ref current
                paramsOptim(i3).Rxn_discharge(i1) = Rxn_discharge_pulse(2); %ref current
      
                %save voltage data
                paramsOptim(i3).mures_charge(i1,:) = interp1((xChargePulse...
                    -xChargePulse(1)),mures_charge_pulse,paramsOptim(i3).t_charge); %ref voltage
                paramsOptim(i3).mures_discharge(i1,:) = interp1((xDischargePulse...
                    -xDischargePulse(1)),mures_discharge_pulse,paramsOptim(i3).t_discharge); %ref voltage
                 paramsOptim(i3).mures_charge_c(i1,:) = interp1((xChargePulse...
                    -xChargePulse(1)),mures_charge_pulse_c,paramsOptim(i3).t_charge); %ref voltage
                 paramsOptim(i3).mures_charge_a(i1,:) = interp1((xChargePulse...
                    -xChargePulse(1)),mures_charge_pulse_a,paramsOptim(i3).t_charge); %ref voltage
                 paramsOptim(i3).mures_discharge_c(i1,:) = interp1((xDischargePulse...
                    -xDischargePulse(1)),mures_discharge_pulse_c,paramsOptim(i3).t_discharge); %ref voltage   
                 paramsOptim(i3).mures_discharge_a(i1,:) = interp1((xDischargePulse...
                    -xDischargePulse(1)),mures_discharge_pulse_a,paramsOptim(i3).t_discharge); %ref voltage   

            end

%             i1
%             i2
%             i3

            % we want to store this data at the t_charge data for the j = 1
            % values
            %storing current values in the dataset
            % ycharge_pulse is voltage
            % Rxn_charge_pulse is the current applied\
            dataRef(i3).y_charge_pulse(i1,i2,:) = interp1((xChargePulse...
                -xChargePulse(1)),mures_charge_pulse,paramsOptim(i3).t_charge); %current         
            dataRef(i3).y_discharge_pulse(i1,i2,:) = interp1((xDischargePulse...
                -xDischargePulse(1)),mures_discharge_pulse,paramsOptim(i3).t_discharge); %current
            
            %storing the voltage at t=0 of relaxation
            dataRef(i3).charge_rlx(i1,i2) = Rxn_charge_pulse(2); %this is voltage difference for relax
            %this value should be negative because we are relaxing
            dataRef(i3).discharge_rlx(i1,i2) = Rxn_discharge_pulse(2);
        
        
        
         end
    end
    
end
end





function [zPulse, zRelax] = invertHPPC()
%perform the optimization

%uncomment this to load single particle simulation data [paramsOptim,
%dataRef, params_nondim, params] = loadDataSingleParticle() uncomment this
%to load MPET simulations data. features many particles as
% %well as a particle size distribution 
[paramsOptim, dataRef, params_nondim, params] = loadDataMPETSimsFullCell()



% define small parameter
eps0 = sqrt(eps);

%choose the unscaled or rescaled norm
% params_nondim.optimize = "together_half_cell"; %separately or together
params_nondim.optimize = "together_full_cell"; %separately or together
params_nondim.norm = "rescaled"; %unscaled or rescaled
params_nondim.norm_type = "L2"; % L2_one_point or L2_all_time

%used for full cell: exact or linear approximation
params_nondim.approx = "exact";
%params_nondim.approx = "linear";


switch params_nondim.optimize 
%     case "separately"
%         % option 1: optimize each SOC separately optimize over each state
%         % of charge z represents the three parameters that are being
%         % optimized for [R_f, c_tilde, c_lyte] initial conditions are [0,
%         % 1, 1]
%         
%         % N is number of cycles for each simulation.
%         N = params.num_cycles;
%         zPulse = zeros(N, 10, 3);
%         zRelax = ones(N, 10);
%         for i1 = 1:N
%             zPulse(i1,1,:) = [eps0,1-eps0,1-eps0];
%         end
% 
%         % do optimization
%         for i1 = 1:N
%             for i2 = 2:10
%                 zPulse(i1, i2, :) = optimizePulse(zPulse(i1, i2-1, :), ...
%                     paramsOptim(i1), i2, ...
%                     params_nondim, dataRef(i1));
%                 zRelax(i1,i2) = optimizeRelax(zRelax(i2-1), ...
%                     paramsOptim(i1), i2, params_nondim, dataRef(i1));
% 
%             end
%         end
% 
%         xMu = squeeze(mean(zPulse));
%         xVar = squeeze(var(zPulse));
% 
%         %fig 1: R_f, fig2: c_tilde, fig3: c_lyte
%         figure(1)
%         for i1 = 1:3
% 
%             subplot(2,2,i1)
%             errorbar(xMu(:,i1), xVar(:,i1), 'DisplayName', 'Fitted Solution')
%             hold on
%             xlabel('Degradation State')
%             if i1 == 1
%                 plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
%                 ylabel('R_f (nondim)')
%             elseif i1 == 2
%                 plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
%                 ylabel('c_{tilde}')
%             elseif i1 == 3
%                 plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
%                 ylabel('c_{+} (M)')
%             end
%             legend()
%         end
%         saveas(gcf,'pulse_fit_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type +'.png')
% 
% 
%         figure(2)
%         plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution')
%         hold on
%         errorbar(mean(zRelax), var(zRelax), 'DisplayName', 'Fitted Solution')
%         xlabel('State')
%         ylabel('c_{+} (M)')
%         legend()
%         saveas(gcf,'relax_fit_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type +'.png')
% 
    case "together"

        % option 2: optimize over all SOC at the same time
        %uncomment the objective function integratedAreaAllSOC
        % z represents the three parameters that are being optimized for
        % [R_f, c_tilde, c_lyte] initial conditions are [0, 1, 1]
        zPulse = zeros(10, 3);
        zRelax = ones(10, 1);
        zPulse(:,1) = eps0*ones(10,1);
        zPulse(:,2) = (1-eps0)*ones(10,1);
        zPulse(:,3) = (1-eps0)*ones(10,1);

        % do optimization
        for i2 = 2:10
            i2
            zPulse(i2, :) = optimizePulse(zPulse(i2-1, :), ...
                paramsOptim, i2, ...
                params_nondim, dataRef);
%             zRelax(i2)   = optimizeRelax(zRelax(i2-1), ...
%                 paramsOptim, i2, params_nondim, dataRef);
        end


        %fig 1: R_f, fig2: c_tilde, fig3: c_lyte
        figure(1)
        for i1 = 1:3

            subplot(2,2,i1)
            plot(zPulse(:,i1),'DisplayName', 'Fitted Solution')
            hold on
            xlabel('Degradation State')
            if i1 == 1
                plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('R_f (nondim)')
            elseif i1 == 2
                plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('c_{tilde}')
            elseif i1 == 3
                plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('c_{+} (M)')
            end
            legend()
        end
        saveas(gcf,'pulse_fit_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type +'.png')

        figure(2)
        plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution')
        hold on
        plot(zRelax, 'DisplayName', 'Fitted Solution')
        xlabel('State')
        ylabel('c_{+} (M)')
        legend()
        saveas(gcf,'relax_fit_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type +'.png')
        save('pulse_fit_bath_R_f_only_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type + '_' + params_nondim.approx + '.mat')

    case "together_full_cell"

        % option 2: optimize over all SOC at the same time
        %uncomment the objective function integratedAreaAllSOC
        % z represents the three parameters that are being optimized for
        % [R_f, c_tilde, c_lyte] initial conditions are [0, 1, 1]
        zPulse = zeros(10, 5);
        zRelax = ones(10, 1);
        zPulse(:,1) = eps0*ones(10,1);
        zPulse(:,2) = (1-eps0)*ones(10,1);
        zPulse(:,3) = (1-eps0)*ones(10,1);
        zPulse(:,4) = eps0*ones(10,1);
        zPulse(:,5) = (1-eps0)*ones(10,1);
        saved_times = zeros(10-1-1,1);

        % do optimization
        for i2 = 2:10
            i2 
            tic()
            zPulse(i2, :) = optimizePulseFullCell(zPulse(i2-1, :), ...
                paramsOptim, i2, ...
                params_nondim, dataRef);         
            saved_times(i2-1) = toc();
%             zRelax(i2)   = optimizeRelaxFullCell(zRelax(i2-1), ...
%                 paramsOptim, i2, params_nondim, dataRef);
        end

        disp("Averaged Time: " + string(mean(saved_times)) + "s")


        %fig 1: R_f, fig2: c_tilde, fig3: c_lyte
        figure(1)
        for i1 = 1:5

            subplot(3,2,i1)
            plot(zPulse(:,i1),'DisplayName', 'Fitted Solution')
            hold on
            xlabel('Degradation State')
            if i1 == 1
                plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('R_{f,c} (nondim)')
            elseif i1 == 2
                plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('c_{tilde,c}')
            elseif i1 == 3
                plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('c_{+} (M)')
            elseif i1 == 4
                plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('R_{f,a} (nondim)')
            elseif i1 == 5
                plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                ylabel('c_{tilde,a}')
               
            end
            legend()
%         saveas(gcf,'full_cell_pulse_fit_' + params_nondim.norm + '_' + params_nondim.optimize + '_' + params_nondim.norm_type +'.png')
        save('current_pulse_results/full_cell_pulse_fit_c_lyte_only_' + params_nondim.norm + '_' + ...
            params_nondim.optimize + '_' + params_nondim.norm_type + '_' +...
            params_nondim.approx + '.mat')


        end

   case "together_half_cell"

        % option 2: optimize over all SOC at the same time
        %uncomment the objective function integratedAreaAllSOC
        % z represents the three parameters that are being optimized for
        % [R_f, c_tilde, c_lyte] initial conditions are [0, 1, 1]
        zPulse = zeros(10, 3);
        zRelax = ones(10, 1);
        zPulse(:,1) = eps0*ones(10,1);
        zPulse(:,2) = (1-eps0)*ones(10,1);
        zPulse(:,3) = (1-eps0)*ones(10,1);
        saved_times = zeros(10-1-1,1);

        % do optimization
        for i2 = 2:8
            i2 
            tic()
            zPulse(i2, :) = optimizePulseHalfCell(zPulse(i2-1, :), ...
                paramsOptim, i2, ...
                params_nondim, dataRef);         
            saved_times(i2-1) = toc();
%             zRelax(i2)   = optimizeRelaxFullCell(zRelax(i2-1), ...
%                 paramsOptim, i2, params_nondim, dataRef);
        end

        disp("Averaged Time: " + string(mean(saved_times)) + "s")

        save('half_cell_pulse_fit_all_deg_' + params_nondim.norm + '_' + ...
            params_nondim.optimize + '_' + params_nondim.norm_type + '_' +...
            params_nondim.approx + '.mat')
        %fig 1: R_f, fig2: c_tilde, fig3: c_lyte
         figure(1)
         for i1 = 1:5
 
             subplot(3,2,i1)
             plot(zPulse(:,i1),'DisplayName', 'Fitted Solution')
             hold on
             xlabel('Degradation State')
             if i1 == 1
                 plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
                 ylabel('R_{f,c} (nondim)')
             elseif i1 == 2
                 plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                 ylabel('c_{tilde,c}')
             elseif i1 == 3
                 plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                 ylabel('c_{+} (M)')
             elseif i1 == 4
                 plot(linspace(0,0.025,10), 'k--', 'DisplayName', 'Reference Solution');
                 ylabel('R_{f,a} (nondim)')
             elseif i1 == 5
                 plot(linspace(1,0.9,10), 'k--', 'DisplayName', 'Reference Solution');
                 ylabel('c_{tilde,a}')
                
             end
             legend()
         end

end
end


function out = dhelper_fundetaf(eta_f, lmbda)
%derivative of helper function for CIET
out = (eta_f.*exp(-(lmbda - (eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)).^2./(4*lmbda))) ...
    ./((exp(-eta_f) + 1).*(eta_f.^2 + lmbda^(1/2) + 1).^(1/2)) - ...
    (lmbda^(1/2).*pi.^(1/2)*exp(-eta_f).*(erf((lmbda - (eta_f.^2 + ...
    lmbda^(1/2) + 1).^(1/2))./(2*lmbda^(1/2))) - 1))./(exp(-eta_f) + 1).^2;
end



function out = h_fun(c, mures, params_nondim)
%convert rxn*R_f to kBT/e
eta = mu_c(c, params_nondim)-mures;
% eta_f = eta - log(c./c_lyte) + rxn*R_f; i_red = helper_fun(-eta_f,
% lambda); i_ox = helper_fun(eta_f, lambda); out =
% k/sqrt(4*pi*lambda)*(v-c).*(c_lyte.*i_red - c.*i_ox) - rxn;
eta_f = eta + log(1./c);
ecd_extras = 1./sqrt(4.0*pi*params_nondim.lambda);
krd = helper_fun(-eta_f, params_nondim.lambda);
kox = helper_fun(eta_f, params_nondim.lambda);
out = ecd_extras*(krd - kox*c);
end



function out = dideta(c, mures, k0, params)
%dideta function for CIET
muh = mu_c(c, params);
eta = muh - mures;
% etaf = eta - log(c);
etaf = eta - log(c);
lmbda = params.lambda;
out = k0.*(1-c)/sqrt(4*pi*lmbda).*...
     (-dhelper_fundetaf(-etaf, lmbda) - c.*dhelper_fundetaf(etaf, lmbda));
end

function rxn = iredoi(c, mures, params_nondim)
%ired/i fraction used in electrolyte fitness value R is dc/dt, from the
%population balance equation. can be BV or linear
% %set symmetry coefs
muh = mu_c(c, params_nondim);
%for simplicity in calculation
etaf = muh-mures- log(c);
rxn = helper_fun(-etaf, params_nondim.lambda)./(helper_fun(-etaf, params_nondim.lambda)-c.*helper_fun(etaf, params_nondim.lambda));
end

function rxn = ired(c, mures, k0, params_nondim)
%ired/ value used in electrolyte concnetration R is dc/dt, from the
%population balance equation. can be BV or linear set symmetry coefs
muh = mu_c(c, params_nondim);
%for simplicity in calculation
etaf = muh-mures- log(c);
rxn = k0/sqrt(4*pi*params_nondim.lambda)*(1-c).*helper_fun(-etaf, params_nondim ...
    .lambda);
end 

function out = W(c, R_f, c_tilde, c_lyte, mures, params, k0_input, lambda_input)
% input parameters for fitness function: c: solid concentration; R_f:
% (nondim) film resistance; c_tilde: capacity loss scaled relative to
% original value; c_lyte: electrolyte concentration in (M), mures:
% electrolyte potential in kBT; params: param functions; k0: optional input
% parameter for exchange current density. otherwise we use the default one
if params.electrolyte == "SM"
    k0 = 10;
else
    k0 = params.k;
end
%overload exchange current density for a full cell
if exist('k0_input','var')
    k0 = k0_input;
end
if exist('lambda_input','var')
    params.lambda = lambda_input;
end
%fitness value for CIET we need iredoi to be slighlty posivie/nevative to
%obtain the values if we use a stefan-maxwell concentrated electrolyte
%model
di = dideta(c, mures, k0, params);
out = 1./(1-R_f.*di/0.0257).* ...
     (c_tilde-c)./(1-c).*(1-iredoi(c, mures, params).*(1-c_lyte).*therm_fac(c_lyte, params));
 %uncommnet this for single particle
end



function out = U(c, R_f, c_tilde, c_lyte, curr, mures, params, k0_input, lambda_input)
% input parameters for fitness function: c: solid concentration; R_f:
% (nondim) film resistance; c_tilde: capacity loss scaled relative to
% original value; c_lyte: electrolyte concentration in (M), mures:
% electrolyte potential in kBT; params: param functions; k0: optional input
% parameter for exchange current density. otherwise we use the default one
if params.electrolyte == "SM"
    k0 = 10;
else
    k0 = params.k;
end
%overload exchange current density for a full cell
if exist('k0_input','var')
    k0 = k0_input;
end
if exist('lambda_input','var')
    params.lambda = lambda_input;
end
%fitness value for CIET we need iredoi to be slighlty posivie/nevative to
%obtain the values if we use a stefan-maxwell concentrated electrolyte
%model
di = dideta(c, mures, k0, params);
out =(1+curr./mures.*R_f./0.0257).* ...
     (1+di.^-1 .* mures.^-1 .* k0 .* h_fun(c, mures, params) .*(c_tilde-1)) ...
     .*(1+di.^-1 .* (mures-log(c./c_tilde)).^-1 .* k0.*(1-c)./sqrt(4.0*pi*params.lambda).*helper_fun(-(mures-log(c./c_tilde)), params.lambda).*(c_lyte-1));
 %uncommnet this for single particle
end



function out = W_ref_anode(c_lyte, mures_a, params)
% input parameters for fitness function: c: solid concentration; R_f:
% (nondim) film resistance; c_tilde: capacity loss scaled relative to
% original value; c_lyte: electrolyte concentration in (M), mures:
% electrolyte potential in kBT; params: param functions; k0: optional input
% parameter for exchange current density. otherwise we use the default one
out = 1-0.5*(1-coth(mures_a/2))*therm_fac(c_lyte, params).*(1-c_lyte);
% out = 1-0.5*therm_fac(c_lyte, params).*(1-c_lyte);
 %uncommnet this for single particle
end


function out = dideta_ratio(c_c, c_a, k0_c, k0_a, lambda_c, lambda_a, mures_c, mures_a, params)
%gets di_c/di_a for a full cell electrode.
params.electrode = "c"; params.lambda = lambda_c; di_c = dideta(c_c, mures_c, k0_c, params);
params.electrode = "a"; params.lambda = lambda_a; di_a = dideta(c_a, mures_a, k0_a, params);
out = di_c./di_a;
end

function out = dideta_ratio_half_cell(c_c, k0_c, lambda_c, mures_c, mures_a, params)
%gets di_c/di_a for the half cell. thus, the dideta_a is for the anode only
%for a simple BV reaction model.
params.electrode = "c"; params.lambda = lambda_c; di_c = dideta(c_c, mures_c, k0_c, params);
di_a = cosh(mures_a/2);
out = di_c./di_a;
end


function out = therm_fac(c, params)
% valoen reimers thermodynamic factor. c is electrolyte concentration in M
out = 1; %default is dilute model
%thermfac = dln(a_l)/dln(c_l)-> a_l^-1-a_l(1)^-1 = \int_1^c_l c_l^-1 *
%thermfac(c_l) dc_l a_l(c_l) = (\int_1^c_l c^-1 * thermfac(c) + 1 )^(-1)

if isfield(params, "electrolyte")
    if params.electrolyte == "SM"
        tmp = 0.601 - 0.24*c.^(0.5) + 0.982*(1 - 0.0052*(298 - 294))*c.^(1.5);
        out = tmp/(1-0.38);
    end
end
end 


function out = a_l(c, params)
% valoen reimers thermodynamic factor. c is electrolyte concentration in M
out = c; %default is dilute model
%thermfac = ln(a+) = int_1^c TF(c_+)/c_+ dc
if isfield(params, "electrolyte")
    if params.electrolyte == "SM"
        out = exp((601.*log(c.^(1/2)))./310 - (24.*c.^(1/2))/31 + (100164.*c.^(3/2))./96875 - 0.2598);
    end
end
end 


function out = dmudc(c, params)
%here we have the scaled chemical potential summed from the two different
%phases dmudc OCV from Colclasure 2020 NMC532
OCV = dColclasure_OCVdx(c);
%the OCV is in units of V, convert to kBT
out = -params.e/params.kBT*OCV;
end


function OCV = dColclasure_OCVdx(x)
%colclasure OCV function derivative
OCV = -3.640117692001490E+03*14*x.^13.0 + 1.317657544484270E+04*13*x.^12.0...
- 1.455742062291360E+04*12*x.^11.0 - 1.571094264365090E+03*11*x.^10.0...
+ 1.265630978512400E+04*10*x.^9.0 - 2.057808873526350E+03*9*x.^8.0...
- 1.074374333186190E+04*8*x.^7.0 + 8.698112755348720E+03*7*x.^6.0...
- 8.297904604107030E+02*6*x.^5.0 - 2.073765547574810E+03*5*x.^4.0...
+ 1.190223421193310E+03*4*x.^3.0 - 2.724851668445780E+02*3*x.^2.0...
+ 2.723409218042130E+01*2*x - 4.158276603609060E+00 +...
-5.573191762723310E-04*exp(6.560240842659690E+00*x.^4.148209275061330E+01)*6.560240842659690E+00*4.148209275061330E+01.*x.^3.148209275061330E+01;
end


function out = mu_c(c, params_nondim)
%here we have the scaled chemical potential summed from the two different
%phases OCV from Colclasure 2020 NMC532
OCV = Colclasure_OCV(c);
%the OCV is in units of V, convert to kBT
out = -params_nondim.e/params_nondim.kBT*OCV + params_nondim.muR_ref_c;
if params_nondim.electrode == "a"
    %if the input is anode
    muRtheta = -1/params_nondim.kToe*0.12;
    muRhomog = graphite_1param_homog_3(c, 3.4, 1.4);
    out = muRhomog + muRtheta + params_nondim.muR_ref_a;
end
end


function muR = graphite_1param_homog_3(y, Omga, Omgb)
%Helper function with low hysteresis and soft tail
width = 5e-2;
tailScl = 5e-2;
muLtail = -tailScl*1./(y.^(0.85));
muRtail = tailScl*1./((1-y).^(0.85));
muRtail = 1.0e1*step_up(y, 1.0, 0.045);
muLlin = (0.15*Omga*12*(0.40-y.^0.98) ...
          .* step_down(y, 0.49, 0.9*width).*step_up(y, 0.35, width));
muRlin = (0.1*Omga*4*(0.74-y) + 0.90*Omgb).*step_up(y, 0.5, 0.4*width);
muLMod = (0. ...
          + 40*(-exp(-y/0.015)) ...
          + 0.75*(tanh((y-0.17)/0.02) - 1) ...
          + 1.0*(tanh((y-0.22)/0.040) - 1) ...
          ).*step_down(y, 0.35, width);
muR = 0.18 + muLMod + muLtail + muRtail + muLlin + muRlin;
end


function out = step_down(x, xc, delta)
out = 0.5.*(-tanh((x - xc)./delta) + 1);
end


function out = step_up(x, xc, delta)
out = 0.5.*(tanh((x - xc)./delta) + 1);
end


function OCV = Colclasure_OCV(x)
%colclasure OCV function
OCV = (5.314735633000300E+00 +...
-3.640117692001490E+03*x.^14.0 + 1.317657544484270E+04*x.^13.0...
- 1.455742062291360E+04*x.^12.0 - 1.571094264365090E+03*x.^11.0...
+ 1.265630978512400E+04*x.^10.0 - 2.057808873526350E+03*x.^9.0...
- 1.074374333186190E+04*x.^8.0 + 8.698112755348720E+03*x.^7.0...
- 8.297904604107030E+02*x.^6.0 - 2.073765547574810E+03*x.^5.0...
+ 1.190223421193310E+03*x.^4.0 - 2.724851668445780E+02*x.^3.0...
+ 2.723409218042130E+01*x.^2.0 - 4.158276603609060E+00*x +...
-5.573191762723310E-04*exp(6.560240842659690E+00*x.^4.148209275061330E+01));
end


function out = R_value(rxn, c, k, lambda, v, R_f, c_lyte, muh, mures, params_nondim)
eta = muh-mures+rxn*R_f/0.0257;
al = a_l(c_lyte, params_nondim);
eta_f = eta + log(al./c);
ecd_extras = (v-c)./sqrt(4.0*pi*lambda);
krd = k*helper_fun(-eta_f, lambda);
kox = k*helper_fun(eta_f, lambda);
out = ecd_extras*(krd*al - kox*c)-rxn;
end


function out = R(c, k, lambda, v, R_f, c_lyte, muh, mures, params_nondim)
%reaction rate
out = fzero(@(rxn) R_value(rxn, c, k, lambda, v, R_f, c_lyte, muh, mures, params_nondim), 1e-1);
% out = fzero(@(rxn) i_value(rxn, c, k, v, R_f, c_lyte, mures, params_nondim), 1);
end



function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end



function out = solve_i_full(mu_values, c_charge_v_c, R_f_c, c_tilde_c, c_lyte, ...
            params_nondim, k0_c, lambda_c, c_charge_v_a, R_f_a, c_tilde_a, ...
            k0_a, lambda_a, curr_constraint)
%takes input values and calculates exact value for deltaV for a full cell
mures_charge_vd_c = mu_values(1);
mures_charge_vd_a = mu_values(2);
%solves for phi_cathode, phi_anode
params_nondim.electrode = "a"; mua = mu_c(c_charge_v_a, params_nondim);
i_a = R(c_charge_v_a, k0_a, lambda_a, c_tilde_a, R_f_a, c_lyte, mua, mures_charge_vd_a, params_nondim);
params_nondim.electrode = "c"; muc = mu_c(c_charge_v_c, params_nondim);
i_c = R(c_charge_v_c, k0_c, lambda_c, c_tilde_c, R_f_c, c_lyte, muc, mures_charge_vd_c, params_nondim);
% current equality
out(1) = i_a + params_nondim.ca_ratio*i_c;
out(2) = i_c*params_nondim.rescale_c - curr_constraint;
end




function [res, fig] = integratedAreaAllSOCFullCell(z, d, paramsOptim, params_nondim, ...
    dataRef, flagPlot)
% d is the degraded state of charge freq is 1/(t-t0) x is the integrated
% area freq = params.freq;

% sum over all different pulse sizes
res = 0;
V = size(paramsOptim(1).voltage_pulse_range, 2);

R_f_c = z(1); c_tilde_c = z(2); 
c_lyte = z(3);
R_f_a = z(4); c_tilde_a = z(5); 
%for a full cell, we need to fit each parameter for each electrode, as well
%as the overall voltage shift. However, we gain two fitness values to
%optimize over as well.

if flagPlot
    % Plot Testing
    fig.fi = figure('Name', 'Visual Comparison Charge/Discharge', ...
        'Position', [50, 50, 560*2, 420]); clf(fig.fi);
    fig.ax = cell(V,2); fig.pl = cell(V,2,2); fig.lg = cell(V,1,2);
    counterPlot = 0;

end

%1. get the c_charge values ofr cathode, anode; 2. get the k0 values for
%cathode, anode 3. optimize

% %nondim by e*len*c_s_ref, but not t_ref
k0_c = 10; k0_a = 0.2;
lambda_c = 3.78; lambda_a = 5;

%MANUALLY OVERWRITE for testing 
% c_tilde_a = 1; R_f_a = 0; R_f_c = 0;
% c_tilde_c = 1;

options = optimset('Display','off') ;
for i1 = 1:V
    for i3 = 1:6
        %% Compare Charge Currents
        % pull system data
        rxn_charge = squeeze(paramsOptim(i3).Rxn_charge(i1));
        c_charge_v_c = paramsOptim(i3).c_charge_c(i1, 1);
        c_charge_v_a = paramsOptim(i3).c_charge_a(i1, 1);
        mures_charge_vd_c = squeeze(paramsOptim(i3).mures_charge_c(i1,:));
        mures_charge_vd_a = squeeze(paramsOptim(i3).mures_charge_a(i1,:));

        Rxn_charge_pulse = squeeze(dataRef(i3).charge_rlx(i1,:));
        Rxn_discharge_pulse = squeeze(dataRef(i3).discharge_rlx(i1,:));

        if params_nondim.approx == "linear"

            params_nondim.electrode = "c";
            U_charge_test_c = U(c_charge_v_c, R_f_c, c_tilde_c, c_lyte, ...
              Rxn_charge_pulse(1)/params_nondim.rescale_c, mures_charge_vd_c, params_nondim);
            params_nondim.electrode = "a";
            U_charge_test_a = U(c_charge_v_a, R_f_a, c_tilde_a, c_lyte, ...
              -Rxn_charge_pulse(1)/params_nondim.rescale_a, mures_charge_vd_a, params_nondim);
            charge_test = U_charge_test_c.*mures_charge_vd_c - U_charge_test_a.*mures_charge_vd_a;

%         %here for testing that the Delta Phi values are equal 
%         dideta_a = dideta(c_charge_v_a, mures_charge_vd_a, k0_a, params_nondim);
%         deltaV = rxn_charge(1).*(W_charge_c_bar-1)...
%             /(params_nondim.rescale_a*dideta_a+W_ca*params_nondim.rescale_c*diratio*dideta_a);
        elseif params_nondim.approx == "exact"
        % pull/calculate test and true charge rate
            optimum_mures = fsolve(@(delta_V) solve_i_full(delta_V, c_charge_v_c, ...
                R_f_c, c_tilde_c, c_lyte, params_nondim, 10, 3.78, ...
                c_charge_v_a, R_f_a, c_tilde_a, 0.2, 5, rxn_charge), [paramsOptim(i3).mures_charge_c(i1,1), paramsOptim(i3).mures_charge_a(i1,1)], options);
            delta_V_solve = optimum_mures(1)-optimum_mures(2);

            charge_test = delta_V_solve;
        end
        charge_true = squeeze(dataRef(i3).y_charge_pulse(i1,d,:));
        
        % calculate the absolute value of the difference in the test and
        % truth
        switch params_nondim.norm
            case "rescaled"
                charge_diff = abs(charge_test - charge_true)./abs(charge_true);
            case "unscaled"
                charge_diff = abs(charge_test - charge_true);
        end
        charge_diff = charge_diff(~isnan(charge_diff));
        
        %% Transform Charge Time to Frequency
        w_ch = paramsOptim(i3).t_charge;
        % FIRST VALUE OF W IS SUPER BIG BECAUSE IT IS JUST A 1/SMALLVALUE.
        % I have just removed it.
        w_ch = w_ch(~isnan(charge_diff));
        
        %% Compare Discharge Currents
        % pull system data
        rxn_discharge = squeeze(paramsOptim(i3).Rxn_discharge(i1));
        c_discharge_v_c = paramsOptim(i3).c_discharge_c(i1, 1);
        c_discharge_v_a = paramsOptim(i3).c_discharge_a(i1, 1);
        mures_discharge_vd_c = squeeze(paramsOptim(i3).mures_discharge_c(i1));
        mures_discharge_vd_a = squeeze(paramsOptim(i3).mures_discharge_a(i1));

%        % option 2: solve for i_a = i_c for delta mures. this options is a
%       little bit more accurate for dilute solution
        %       little bit more accurate for dilute solution   
%     
%        %test the voltage values
%        dideta_a = dideta(c_discharge_v_a, mures_discharge_vd_a, k0_a, params_nondim);
%        deltaV = rxn_discharge(1).*(W_discharge_c_bar-1)...
%             /(params_nondim.rescale_a*dideta_a+W_ca*params_nondim.rescale_c*diratio*dideta_a);

        % pull/calculate test and true charge rate
        if params_nondim.approx == "linear"

            params_nondim.electrode = "c";
            U_discharge_test_c = U(c_discharge_v_c, R_f_c, c_tilde_c, c_lyte, ...
              Rxn_discharge_pulse(1)/params_nondim.rescale_c, mures_discharge_vd_c, params_nondim);
            params_nondim.electrode = "a";
            U_discharge_test_a = U(c_discharge_v_a, R_f_a, c_tilde_a, c_lyte, ...
              -Rxn_discharge_pulse(1)/params_nondim.rescale_a, mures_discharge_vd_a, params_nondim);
            discharge_test = U_discharge_test_c*mures_discharge_vd_c - U_discharge_test_a*mures_discharge_vd_a;

        elseif params_nondim.approx == "exact"
            optimum_mures = fsolve(@(delta_V) solve_i_full(delta_V, c_discharge_v_c, ...
                R_f_c, c_tilde_c, c_lyte, params_nondim, 10, 3.78, ...
                c_discharge_v_a, R_f_a, c_tilde_a, 0.2, 5, rxn_discharge), [paramsOptim(i3).mures_discharge_c(i1,1), paramsOptim(i3).mures_discharge_a(i1,1)], options);
            delta_V_solve = optimum_mures(1)-optimum_mures(2);


            discharge_test = delta_V_solve;
        end
        discharge_true = squeeze(dataRef(i3).y_discharge_pulse(i1,d,:));
        
        % calculate the absolute value of the difference in the test and
        % truth
        switch params_nondim.norm
            case "rescaled"
                discharge_diff = abs(discharge_test - discharge_true)./abs(discharge_true);
            case "unscaled"
                discharge_diff = abs(discharge_test - discharge_true);
        end
        discharge_diff = discharge_diff(~isnan(discharge_diff));
        
        %% Transform Charge Time to Frequency
        w_di = paramsOptim(i3).t_discharge;
        % FIRST VALUE OF W IS SUPER BIG BECAUSE IT IS JUST A 1/SMALLVALUE.
        % I have just removed it.
        w_di = w_di(~isnan(discharge_diff));
        
        switch params_nondim.norm_type 
            case "L1"
%                 res = res + ...
%                     abs(trapz(w_ch, charge_diff)) + ...
%                     abs(trapz(w_di, discharge_diff));  
                res = res + ...
                    abs(charge_diff(1)) + ...
                    abs(discharge_diff(1));  
            case "L1_charge"
                res = res + ...
                    abs(charge_diff(1));    
            case "L1_discharge"
                res = res + ...
                    abs(discharge_diff(1));    
            case "L2"
                res = res + ...
                    charge_diff(1).^2 + ...
                    discharge_diff(1).^2;  

        end

        
    end
    fig = [];   
   
end
    
end




function [z] = optimizePulseFullCell(z0, paramsOptim, d, ...
    params_nondim, dataRef)
% [Z] = OPTIMIZEPULSE(X0,OPTIMPARAMS,D,PARAMS_NONDIM,DATAREF) optimizes the
% values of (R_f, c_tilde, c_lyte) for a given pulsing strategy.
% 
% Author: Debbie Zhuang Date: 12/05/2022
%
%% Input:
% (1) z0: size: 1x3
%         class: array
%     description: initial guess for the value of the degradation
%                  parameters
% (2) optimParams: size: 1x1
%                  class: struct
%     description: optimization parameters for the system.
% (3) d: size: 1x1
%        class: int
%     description: index of degradation, varies from 1 to 10.
% (4) params_nondim: size: 1x1
%                    class: struct
%     description: contains non-dimensionalized parameters.
% (5) dataRef: size: 1x1
%              class: struct
%     description: contains the reference system data
%
%% Output:
% (1) z: size: 1x3
%        class: array
%     description: optimized degradation parameters
%
%% Do Optimization

% define inequality constraints range: R_f_c >=0, 0.5 <= ctilde <= 1, 0.5 <=
% ctilde <= 1, R_f >=0, 0 <= c_lyte <=1
A = -eye(5); A(6,2) = 1; A(7,3) = 1; A(8,5) = 1; A(9,1) = 1; A(10,4) = 1;
b = [0; -0.5; -0.55; 0; -0.5; 1; 1; 1; 10; 10];

% define optimization options
options = optimset('Display','off') ;

V = size(paramsOptim(1).voltage_pulse_range, 2);

V_ref_charge = zeros(V, 6);
V_ref_discharge = zeros(V, 6);

for i1 = 1:V
    for i3 = 1:6
        %% Compare Charge Currents
        % pull system data
        rxn_charge = squeeze(paramsOptim(i3).Rxn_charge(i1));
        c_charge_v_c = paramsOptim(i3).c_charge_c(i1, 1);
        c_charge_v_a = paramsOptim(i3).c_charge_a(i1, 1);
        mures_charge_vd = squeeze(paramsOptim(i3).mures_charge(i1,:));
        mures_charge_vd_c = squeeze(paramsOptim(i3).mures_charge_c(i1,:));
        mures_charge_vd_a = squeeze(paramsOptim(i3).mures_charge_a(i1,:));
        rxn_discharge = squeeze(paramsOptim(i3).Rxn_discharge(i1));
        c_discharge_v_c = paramsOptim(i3).c_discharge_c(i1, 1);
        c_discharge_v_a = paramsOptim(i3).c_discharge_a(i1, 1);
        mures_discharge_vd = squeeze(paramsOptim(i3).mures_discharge(i1,:));
        mures_discharge_vd_c = squeeze(paramsOptim(i3).mures_discharge_c(i1,:));
        mures_discharge_vd_a = squeeze(paramsOptim(i3).mures_discharge_a(i1,:));

% exact solution only:
% solve for the original voltage minus the expected voltage
%need to solve inversely, mu_c - mu_a of real vs invR(R) - invR(R) of
%expected

%         % option 1: linearized solution
%         optimum_mures = fsolve(@(delta_V) solve_i_full(delta_V, c_charge_v_c, ...
%             0, 1, 1, params_nondim, 10, 3.78, ...
%             c_charge_v_a, 0, 1, 0.2, 5, rxn_charge), [paramsOptim(i3).mures_charge_c(1,1), paramsOptim(i3).mures_charge_a(1,1)], options);
%         delta_V_solve = optimum_mures(1)-optimum_mures(2);
%         
%         % minus reference solution
%         V_ref_charge(i1, i3) = delta_V_solve - mures_charge_vd(1);
% 
%             %       little bit more accurate for dilute solution
%         optimum_mures = fsolve(@(delta_V) solve_i_full(delta_V, c_discharge_v_c, ...
%             0, 1, 1, params_nondim, 10, 3.78, ...
%             c_discharge_v_a, 0, 1, 0.2, 5, rxn_discharge), [paramsOptim(i3).mures_discharge_c(1,1), paramsOptim(i3).mures_discharge_a(1,1)], options);
%         delta_V_solve = optimum_mures(1)-optimum_mures(2);
% 
%         V_ref_discharge(i1, i3) = delta_V_solve - mures_discharge_vd(1);


    end
end

% params_nondim.V_ref_charge = V_ref_charge;
% params_nondim.V_ref_discharge = V_ref_discharge;

options = optimoptions('fmincon', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-8, ...
    'Display', 'iter');


% % define objective function
switch params_nondim.optimize
    case "together_full_cell"
        % if we choose to optimize over all at the same time, we redefine
        % objective function
        objFun = @(z) integratedAreaAllSOCFullCell(z, d, paramsOptim, params_nondim, ...
            dataRef, false);
end

% % optimize! [z, ~] = fmincon(objFun, z0, A, b, [],[],[],[],[], options);
% if we optimize over all, we have
[z, ~] = fmincon(objFun, z0, A, b, [],[],[],[],[], options);
switch params_nondim.optimize
    case "together_full_cell"
        [objVal, fig] = integratedAreaAllSOCFullCell(z, d, paramsOptim, params_nondim, ...
            dataRef, false);
end

fprintf('''''''''''''''''''''''''''''''''''''''''''''''''''''''''\n');
fprintf('Inverted Parameters.\n')
disp('z = '); disp(z);
disp('objVal = '); disp(objVal);
fprintf('''''''''''''''''''''''''''''''''''''''''''''''''''''''''\n');

end


function out = nearZero(x, tol)
out = (abs(x) < tol);
end
