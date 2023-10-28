clear;clc;

maxNumCompThreads(1);

% invert HPPC data

% ??????????????? data loaded format: c data at each discretization, for
% 1:200; chemical potential; reaction rate; indicator pulse data: +i, -i is
% indicator relaxation data: +i+N, -i-N is indicator



[xPulse, xRelax] = invertHPPC();


function vq = interp_colors(map_less, len)
vq = zeros(len, 3);
for i = 1:3
    xq = linspace(1, size(map_less,1), len);
    vq(:,i) = interp1(1:size(map_less,1), map_less(:,i), xq);
end
end


function [paramsOptim, dataRef, params_nondim, params] = loadDataMPETSimsFullCell()
%%% load data from single particle simulations
% returns: loaded current and voltage data for current pulses and
% relaxation segments

%R0_bar = 0.0007350654435686459 ish;


% pull one of the datasets...any of them
load('../HPPC_data_5meV_pulse.mat');

% set voltage params of interest
voltRange = [5, 10, 50, 100];
%voltRange = [150];
%voltRange = [100];

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
clyte = 1000:-11.1:900;
Rf = 0:0.0025:0.026;

load('../HPPC_data_5meV_pulse.mat');

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
params_nondim.rescale_a = L_a * (1-poros_a)*P_L_a*mean(3./particle_sizes);
params_nondim.ca_ratio = params_nondim.rescale_c/params_nondim.rescale_a;
params_nondim.electrolyte = 'SM';
%params_nondim.electrolyte = 'dilute';

for i1 = 1:xVar 
    % voltage pulse size
       
    % guess value at each state degradation
    %for i2 = 1:10
    for i2 = 1:8
%degradation state
          load('C:\Users\ChemeGrad2019\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\dezhuang\projects\hppc_data\MPET_full_cell\c_tilde_c_only\voltage_pulses_'+ ...
            string(paramsOptim(1).voltage_pulse_range(i1)*1e-3) + ...
            'V\ctilde_' + sprintf('%.4f',ctilde(i2)) + ...
            '\sim_output\output_data.mat');
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
                mures_charge_pulse_c = -mean(phi_bulk_c(idxCharge,:),2) + mean(phi_lyte_c(idxCharge,:),2) - log(a_l(c_lyte_c(idxCharge), params_nondim));
                mures_charge_pulse_a = -mean(phi_bulk_a(idxCharge,:),2) + mean(phi_lyte_a(idxCharge,:),2) - log(a_l(c_lyte_c(idxCharge), params_nondim));
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
                mures_discharge_pulse_c = -mean(phi_bulk_c(idxDischarge,:),2) + mean(phi_lyte_c(idxDischarge,:),2) - log(a_l(c_lyte_c(idxDischarge), params_nondim));
                mures_discharge_pulse_a = -mean(phi_bulk_a(idxDischarge,:),2) + mean(phi_lyte_a(idxDischarge,:),2) - log(a_l(c_lyte_c(idxDischarge), params_nondim));
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
                paramsOptim(i3).mures_charge = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_discharge = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_charge_c = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_charge_a = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_discharge_c = zeros(length(paramsOptim(1).voltage_pulse_range),1);
                paramsOptim(i3).mures_discharge_a = zeros(length(paramsOptim(1).voltage_pulse_range),1);
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
             
                %store the time data for charge/discharge in the
                %optim_params
               
                paramsOptim(i3).Rxn_charge(i1,:) = interp1((xChargePulse...
                    -xChargePulse(1)),Rxn_charge_pulse,paramsOptim(i3).t_charge); %ref current
                paramsOptim(i3).Rxn_discharge(i1,:) = interp1((xDischargePulse...
                    -xDischargePulse(1)),Rxn_discharge_pulse,paramsOptim(i3).t_discharge); %ref current
      
                %save voltage data
                paramsOptim(i3).mures_charge(i1) = mures_charge_pulse(2); %ref voltage
                paramsOptim(i3).mures_discharge(i1) = mures_discharge_pulse(2); %ref voltage
                paramsOptim(i3).mures_charge_c(i1) = mures_charge_pulse_c(2); %ref voltage
                paramsOptim(i3).mures_charge_a(i1) = mures_charge_pulse_a(2); %ref voltage
                paramsOptim(i3).mures_discharge_c(i1) = mures_discharge_pulse_c(2); %ref voltage   
                paramsOptim(i3).mures_discharge_a(i1) = mures_discharge_pulse_a(2); %ref voltage   

                %after the charge/discharge pulse
                paramsOptim(i3).charge_rlx(i1) = mures_charge_rlx(2); %ref voltage
                paramsOptim(i3).discharge_rlx(i1) = mures_discharge_rlx(2); % ref votlage

            end

       
            % we want to store this data at the t_charge data for the j = 1
            % values
            %storing current values in the dataset
            dataRef(i3).y_charge_pulse(i1,i2,:) = interp1((xChargePulse...
                -xChargePulse(1)),Rxn_charge_pulse,paramsOptim(i3).t_charge); %current
            dataRef(i3).y_discharge_pulse(i1,i2,:) = interp1((xDischargePulse...
                -xDischargePulse(1)),Rxn_discharge_pulse,paramsOptim(i3).t_discharge); %current
            
            %storing the voltage at t=0 of relaxation
            dataRef(i3).charge_rlx(i1,i2) = mures_charge_rlx(2)-mures_charge_pulse(end); %this is voltage difference for relax
            %this value should be negative because we are relaxing
            dataRef(i3).discharge_rlx(i1,i2) = mures_discharge_rlx(2)-mures_discharge_pulse(end);
        
        
         end
    end
    
end
end



function [zPulse, zRelax] = invertHPPC()
%perform the optimization

[paramsOptim, dataRef, params_nondim, params] = loadDataMPETSimsFullCell()




map_pink_less = [252,197,192
250,159,181
247,104,161
221,52,151
174,1,126
122,1,119
73,0,106]/255;
map_pink = interp_colors(map_pink_less, 50);

map_blue_less = [247,252,240
224,243,219
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,104,172
8,64,129]/255;
map_blue = interp_colors(map_blue_less, 50);


map_green_less = [
    204,236,230
153,216,201
102,194,164
65,174,118
35,139,69
0,109,44
0,68,27]/255;
map_green = interp_colors(map_green_less, 50);
map_green = (map_green + map_blue)/2;
%map_green = (map_blue + map_purple)/2;

map_orange_less = [254,227,145
254,196,79
254,153,41
236,112,20
204,76,2
153,52,4
102,37,6]/255;
map_orange = interp_colors(map_orange_less, 50);

map_purple_less = [224,236,244
191,211,230
158,188,218
140,150,198
140,107,177
136,65,157
129,15,124
77,0,75]/255;
map_purple = interp_colors(map_purple_less, 50);




% define small parameter
eps0 = sqrt(eps);

%choose the unscaled or rescaled norm
params_nondim.norm = "rescaled"; %unscaled or rescaled
params_nondim.optimize = "together_full_cell"; %separately or together
params_nondim.norm_type = "L1"; % L1 or L2
% params_nondim.norm = "unscaled"; %unscaled or rescaled
% params_nondim.optimize = "together"; %separately or together or together_full_cell
% params_nondim.norm_type = "L2"; % L1 or L2

%used for full cell: exact or linear approximation
% params_nondim.approx = "exact";
params_nondim.approx = "linear";

switch params_nondim.optimize 
    case "together_full_cell"

        integratedAreaAllSOCFullCell(paramsOptim, params_nondim, map_blue_less, map_green_less,...
            dataRef);       
end
end


function out = dhelper_fundetaf(eta_f, lmbda)
%derivative of helper function for CIET
out = (eta_f.*exp(-(lmbda - (eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)).^2./(4*lmbda))) ...
    ./((exp(-eta_f) + 1).*(eta_f.^2 + lmbda^(1/2) + 1).^(1/2)) - ...
    (lmbda^(1/2).*pi.^(1/2)*exp(-eta_f).*(erf((lmbda - (eta_f.^2 + ...
    lmbda^(1/2) + 1).^(1/2))./(2*lmbda^(1/2))) - 1))./(exp(-eta_f) + 1).^2;
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


function out = dideta_ratio(c_c, c_a, k0_c, k0_a, lambda_c, lambda_a, mures_c, mures_a, params)
%gets di_c/di_a at the not rescale value
params.electrode = "c"; params.lambda = lambda_c; di_c = dideta(c_c, mures_c, k0_c, params);
params.electrode = "a"; params.lambda = lambda_a; di_a = dideta(c_a, mures_a, k0_a, params);
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
%thermfac = dln(a_l)/dln(c_l)-> a_l^-1-a_l(1)^-1 = \int_1^c_l c_l^-1 *
%thermfac(c_l) dc_l a_l(c_l) = (\int_1^c_l c^-1 * thermfac(c) dc)^(-1)
if isfield(params, "electrolyte")
    if params.electrolyte == "SM"
        out =  exp((601.*log(c.^(1/2)))./310 - (24.*c.^(1/2))/31 + (100164.*c.^(3/2))./96875 - 0.2598);
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
% c_eff = c./(1-v); %scale by # of available sites %c_eff = c; eta = (muh -
% mures); alpha = 0.5; % rxn = params_nondim.k.*(165*c_eff.^5 -
% 752.4*c_eff.^4 + 1241*c_eff.^3 ... %       - 941.7*c_eff.^2 + 325*c_eff -
% 35.85).*(exp(-alpha*eta)-exp((1-alpha)*eta)); rxn =
% params_nondim.k.*sqrt(c_eff.*(1-c_eff)).*(exp(-alpha*eta)-exp((1-alpha)*eta));

%convert rxn*R_f to kBT/e
eta = muh-mures+rxn*R_f/0.0257;
% eta_f = eta - log(c./c_lyte) + rxn*R_f; i_red = helper_fun(-eta_f,
% lambda); i_ox = helper_fun(eta_f, lambda); out =
% k/sqrt(4*pi*lambda)*(v-c).*(c_lyte.*i_red - c.*i_ox) - rxn;
al = a_l(c_lyte, params_nondim);
eta_f = eta + log(al./c);
ecd_extras = (v-c)./sqrt(4.0*pi*lambda);
krd = k*helper_fun(-eta_f, lambda);
kox = k*helper_fun(eta_f, lambda);
out = ecd_extras*(krd*al - kox*c)-rxn;
end


function out = R(c, k, lambda, v, R_f, c_lyte, muh, mures, params_nondim)
out = fzero(@(rxn) R_value(rxn, c, k, lambda, v, R_f, c_lyte, muh, mures, params_nondim), 1e-1);
end


function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end




function W_charge = solve_W(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
    c_c, c_a, mures_c, mures_a, k0_c, k0_a, lambda_c, lambda_a)
%takes input values and calculates exact value for deltaV
params_nondim.electrode = "c";
W_charge_c_bar = W(c_c, R_f_c, c_tilde_c, c_lyte, ...
    mures_c, params_nondim, k0_c, lambda_c);
%shoot, we need to get the anode SOC too
params_nondim.electrode = "a";
W_charge_a_bar = W(c_a, R_f_a, c_tilde_a, c_lyte, ...
    mures_a, params_nondim, k0_a, lambda_a);
diratio = dideta_ratio(c_c, c_a, k0_c, k0_a, ...
    lambda_c, lambda_a, mures_c, mures_a, params_nondim);
%             %linearized solution. can also solve directly
 W_charge = (W_charge_c_bar + W_charge_a_bar.*params_nondim.ca_ratio.*diratio)./ ...
     (1+params_nondim.ca_ratio.*diratio);
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
eta = mu_c(c, params)-mures;
%out =(1+curr./eta.*R_f*0.0257).* ...
%charge discharge does not flip curr/eta sign, U should not change signs
%charge discharge flips 1/(di*eta) sign, but h also flips sign so shouldn't
%change signs
%electrolyte should flip signs 
out =(1-curr./eta.*R_f/0.0257).* ...
     (1-di.^-1 .* eta.^-1 .* k0/sqrt(4*pi*params.lambda) .* h_fun(c, mures, params) .*(c_tilde-1)) ...
     .*(1+di.^-1 .* (eta-log(c./c_tilde)).^-1 .* k0.*(1-c)./sqrt(4.0*pi*params.lambda).*helper_fun(-(mures-log(c./c_tilde)), params.lambda).*(c_lyte-1));
% disp("curr value")
% curr
% disp("calc value")
% R(c, k0, params.lambda, 1, 0, 1, mu_c(c, params), mures, params)
% out =(1+curr./eta.*R_f/0.0257).* ...
%      (1+di.^-1 .* eta.^-1 .* curr./(c_tilde-c) .*(c_tilde-1)) ...
%      .*(1+di.^-1 .* (eta-log(c./c_tilde)).^-1 .* k0.*(1-c)./sqrt(4.0*pi*params.lambda).*helper_fun(-(mures-log(c./c_tilde)), params.lambda).*(c_lyte-1));
%uncommnet this for single particle
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



function W_charge = solve_W_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
    c_c, c_a, k0_c, k0_a, lambda_c, lambda_a, mures_c, mures_a, curr)
%takes input values and calculates exact value for deltaV
params_nondim.electrode = "c";
%film ressitance may be wrong, or current may be wrong. fil resistance
%being wrong seems unlikely

U_charge_test_c = U(c_c, R_f_c, c_tilde_c, c_lyte, ...
  curr/params_nondim.rescale_c, mures_c, params_nondim, k0_c, lambda_c)
eta_c = mu_c(c_c, params_nondim)-mures_c;
params_nondim.electrode = "a";
%when current is negative, does U change positive or negatively?
U_charge_test_a = U(c_a, R_f_a, c_tilde_a, c_lyte, ...
  -curr/params_nondim.rescale_a, mures_a, params_nondim, k0_a, lambda_a)
eta_a = mu_c(c_a, params_nondim)-mures_a;
W_charge = U_charge_test_c.*mures_c - U_charge_test_a.*mures_a;

W_charge = W_charge/(mures_c-mures_a);
end



function out = solve_i_full(delta_V, c_charge_v_c, R_f_c, c_tilde_c, c_lyte, ...
            mures_charge_vd_c, params_nondim, k0_c, lambda_c, c_charge_v_a, R_f_a, c_tilde_a, ...
            mures_charge_vd_a, k0_a, lambda_a)
%takes input values and calculates exact value for deltaV
params_nondim.electrode = "a"; mua = mu_c(c_charge_v_a, params_nondim);
i_a = R(c_charge_v_a, k0_a, lambda_a, c_tilde_a, R_f_a, c_lyte, mua, mures_charge_vd_a + delta_V, params_nondim);
params_nondim.electrode = "c"; muc = mu_c(c_charge_v_c, params_nondim);
i_c = R(c_charge_v_c, k0_c, lambda_c, c_tilde_c, R_f_c, c_lyte, muc, mures_charge_vd_c + delta_V, params_nondim);
out = i_a + params_nondim.ca_ratio*i_c;
end



function out = solve_W_exact(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
    c_c, c_a, mures_c, mures_a, k0_c, k0_a, lambda_c, lambda_a)
%takes input values and calculates exact value for deltaV
if length(R_f_c) == 1
    R_f_c = R_f_c*ones(25,1);
end
if length(R_f_a) == 1
    R_f_a = R_f_a*ones(25,1);
end
if length(c_tilde_c) == 1
    c_tilde_c = c_tilde_c*ones(25,1);
end
if length(c_tilde_a) == 1
    c_tilde_a = c_tilde_a*ones(25,1);
end
if length(c_lyte) == 1
    c_lyte = c_lyte*ones(25,1);
end
delta_V = ones(25,1);
out = ones(25,1);

for j = 1:25
    delta_V(j) = fsolve(@(delta_V) solve_i_full(delta_V, c_c, ...
        R_f_c(j), c_tilde_c(j), c_lyte(j), ...
        mures_c + log(a_l(c_lyte(j), params_nondim)), params_nondim, k0_c, lambda_c, ...
        c_a, R_f_a(j), c_tilde_a(j), mures_a + log(a_l(c_lyte(j), params_nondim)), k0_a, ...
        lambda_a), 0);
    params_nondim.electrode = "a"; mua = mu_c(c_a, params_nondim);
    i_a = R(c_a, k0_a, lambda_a, c_tilde_a(j), R_f_a(j), c_lyte(j), mua, mures_a + log(a_l(c_lyte(j), params_nondim)) + delta_V(j), params_nondim);
    i_a_ref = R(c_a, k0_a, lambda_a, 1, 0, 1, mua, mures_a, params_nondim);
    out(j) = i_a/i_a_ref;
end
out = out/out(1); %divide by reference value;
end


function out = solve_W_exact_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
    c_c, c_a, k0_c, k0_a, lambda_c, lambda_a, curr)
%takes input values and calculates exact value for deltaV
if length(R_f_c) == 1
    R_f_c = R_f_c*ones(25,1);
end
if length(R_f_a) == 1
    R_f_a = R_f_a*ones(25,1);
end
if length(c_tilde_c) == 1
    c_tilde_c = c_tilde_c*ones(25,1);
end
if length(c_tilde_a) == 1
    c_tilde_a = c_tilde_a*ones(25,1);
end
if length(c_lyte) == 1
    c_lyte = c_lyte*ones(25,1);
end
V_c = ones(25,1);
V_a = ones(25,1);
V_c_ref = ones(25,1);
V_a_ref = ones(25,1);
out = ones(25,1);

for j = 1:25
    params_nondim.electrode = "a"; muaguess = mu_c(c_a, params_nondim);
    params_nondim.electrode = "c"; mucguess = mu_c(c_c, params_nondim);
    mus = fsolve(@(delta_V) solve_i_full_current(delta_V, c_c, ...
        R_f_c(j), c_tilde_c(j), c_lyte(j), ...
        params_nondim, k0_c, lambda_c, ...
        c_a, R_f_a(j), c_tilde_a(j), k0_a, ...
        lambda_a, curr), [mucguess, muaguess]);
    V_c(j) = mus(1); V_a(j) = mus(2);
    mus = fsolve(@(delta_V) solve_i_full_current(delta_V, c_c, ...
        0, 1, 1, ...
        params_nondim, k0_c, lambda_c, ...
        c_a, 0, 1, k0_a, ...
        lambda_a, curr), [mucguess, muaguess]);
    V_c_ref(j) = mus(1); V_a_ref(j) = mus(2);
    out(j) = (V_c(j)-V_a(j))/(V_c_ref(j)-V_a_ref(j));
end
%out = out/out(1); %divide by reference value;
end



function out = solve_i_full_current(mu_values, c_charge_v_c, R_f_c, c_tilde_c, c_lyte, ...
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





function integratedAreaAllSOCFullCell(paramsOptim, params_nondim, map_blue_less, map_green_less, ...
    dataRef)
% d is the degraded state of charge freq is 1/(t-t0) x is the integrated
% area freq = params.freq;

% sum over all different pulse sizes
% %nondim by e*len*c_s_ref, but not t_ref
k0_c = 10; k0_a = 0.2;
lambda_c = 3.78; lambda_a = 5;

%MANUALLY OVERWRITE for testing 
% c_tilde_a = 1; R_f_a = 0; R_f_c = 0;
% c_tilde_c = 1;

options = optimset('Display','off') ;
i1 = 4; i3 = 3;
%% Compare Charge Currents
% pull system data
rxn_charge = squeeze(paramsOptim(i3).Rxn_charge(i1,:));
c_charge_v_c = paramsOptim(i3).c_charge_c(i1, 1);
c_charge_v_a = paramsOptim(i3).c_charge_a(i1, 1);
mures_charge_vd_c = squeeze(paramsOptim(i3).mures_charge_c(i1));
mures_charge_vd_a = squeeze(paramsOptim(i3).mures_charge_a(i1));

%% Compare Discharge Currents
% pull system data
rxn_discharge = squeeze(paramsOptim(i3).Rxn_discharge(i1, :));
c_discharge_v_c = paramsOptim(i3).c_discharge_c(i1, 1);
c_discharge_v_a = paramsOptim(i3).c_discharge_a(i1, 1);
mures_discharge_vd_c = squeeze(paramsOptim(i3).mures_discharge_c(i1));
mures_discharge_vd_a = squeeze(paramsOptim(i3).mures_discharge_a(i1));


for j = 1:5
    R_f_c = 0; c_tilde_c = 1; c_lyte = 1; R_f_a = 0; c_tilde_a = 1;
    figure(j)
    if j == 1
        R_f_c = linspace(0,0.2,25);
        deg = R_f_c;
    elseif j == 2
        c_tilde_c = linspace(1,0.8,25);
        deg = c_tilde_c;
    elseif j == 3
        c_lyte = linspace(1,0.8,25);
        deg = c_lyte;
    elseif j == 4
        R_f_a = linspace(0,0.2,25);
        deg = R_f_a;
    elseif j == 5
        c_tilde_a = linspace(1,0.8,25);
        deg = c_tilde_a;
    end
   

    %do fro small and large overpotential
    W_charge = solve_W(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, mures_charge_vd_c(1), mures_charge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_discharge = solve_W(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, mures_discharge_vd_c(1), mures_discharge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_charge_exact = solve_W_exact(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, mures_charge_vd_c(1), mures_charge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_discharge_exact = solve_W_exact(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, mures_discharge_vd_c(1), mures_discharge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    voltage_charge_exact = solve_W_exact_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, k0_c, k0_a, lambda_c, lambda_a,1)
    voltage_discharge_exact = solve_W_exact_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, k0_c, k0_a, lambda_c, lambda_a,-1)
    voltage_charge = solve_W_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, k0_c, k0_a, lambda_c, lambda_a, mures_charge_vd_c(1), mures_charge_vd_a(1), rxn_charge(1))
    voltage_discharge = solve_W_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, k0_c, k0_a, lambda_c, lambda_a, mures_discharge_vd_c(1), mures_discharge_vd_a(1), rxn_discharge(1))
    
    if j ~= 1
%        figure(j)
%        plot(deg, W_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'HandleVisibility', 'off')
%        hold on
%        plot(deg, W_discharge, 'LineWidth',2,'Color', map_blue_less(end,:), 'HandleVisibility', 'off')
        figure(j+5)
        plot(deg, W_charge_exact, 'LineWidth',2,'Color', map_green_less(end-3,:))
        hold on
        plot(deg, W_discharge_exact, 'LineWidth',2,'Color', map_green_less(end,:))
        plot(deg, W_charge, 'o', 'LineWidth',2,'Color', map_green_less(end-3,:),'HandleVisibility', 'off')
        plot(deg, W_discharge, 'o', 'LineWidth',2,'Color', map_green_less(end,:), 'HandleVisibility', 'off')
        % plot(deg, voltage_charge_exact, 'LineWidth',2,'Color', map_blue_less(end-3,:))
        % plot(deg, voltage_discharge_exact, 'LineWidth',2,'Color', map_blue_less(end,:))
        % plot(deg, voltage_charge, 'o', 'LineWidth',2,'Color', map_blue_less(end-3,:),'HandleVisibility', 'off')
        % plot(deg, voltage_discharge, 'o', 'LineWidth',2,'Color', map_blue_less(end,:), 'HandleVisibility', 'off')

        if j == 5
            legend('Current Discharge', 'Current Charge','Fontsize', 16, 'Location', 'southwest')
            legend boxoff
        end
        set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
        export_fig(sprintf('fig_sensitivity%d.png', j+5), '-transparent', '-m3', '-r300')

        
    else
 %       figure(j)
 %       plot(deg, W_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:))
 %       hold on
 %       plot(deg, W_discharge, 'LineWidth',2,'Color', map_blue_less(end,:))
        % legend('Charge', 'Discharge', 'Fontsize', 16, 'Location', 'southwest')
        % legend boxoff
        figure(j+5)
        h1 = plot(deg, W_charge_exact, 'LineWidth',2,'Color', map_green_less(end-3,:))
        hold on
        plot(deg, W_discharge_exact, 'LineWidth',2,'Color', map_green_less(end,:),'HandleVisibility', 'off')
        h2 = plot(deg, W_charge, 'o', 'LineWidth',2,'Color', map_green_less(end-3,:))
        plot(deg, W_discharge, 'o', 'LineWidth',2,'Color', map_green_less(end,:), 'HandleVisibility', 'off')
        % plot(deg, voltage_charge_exact, 'LineWidth',2,'Color', map_blue_less(end-3,:))
        % plot(deg, voltage_discharge_exact, 'LineWidth',2,'Color', map_blue_less(end,:))
        % plot(deg, voltage_charge, 'o', 'LineWidth',2,'Color', map_blue_less(end-3,:),'HandleVisibility', 'off')
        % plot(deg, voltage_discharge, 'o', 'LineWidth',2,'Color', map_blue_less(end,:), 'HandleVisibility', 'off')        
        legend([h1, h2], {'Exact', 'Linear'}, 'Fontsize', 16, 'Location', 'southwest')
        legend boxoff
        set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
        export_fig(sprintf('fig_sensitivity%d.png', j+5), '-transparent', '-m3', '-r300')

    end

%    set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
%    export_fig(sprintf('fig_sensitivity%d.png', j), '-transparent', '-m3', '-r300')

        
end


i1 = 4; i3 = 3;
%% Compare Charge Currents
% pull system data
rxn_charge = squeeze(paramsOptim(i3).Rxn_charge(i1,:));
c_charge_v_c = paramsOptim(i3).c_charge_c(i1, 1);
c_charge_v_a = paramsOptim(i3).c_charge_a(i1, 1);
mures_charge_vd_c = squeeze(paramsOptim(i3).mures_charge_c(i1));
mures_charge_vd_a = squeeze(paramsOptim(i3).mures_charge_a(i1));

%% Compare Discharge Currents
% pull system data
rxn_discharge = squeeze(paramsOptim(i3).Rxn_discharge(i1, :));
c_discharge_v_c = paramsOptim(i3).c_discharge_c(i1, 1);
c_discharge_v_a = paramsOptim(i3).c_discharge_a(i1, 1);
mures_discharge_vd_c = squeeze(paramsOptim(i3).mures_discharge_c(i1));
mures_discharge_vd_a = squeeze(paramsOptim(i3).mures_discharge_a(i1));


for j = 1:5
    R_f_c = 0; c_tilde_c = 1; c_lyte = 1; R_f_a = 0; c_tilde_a = 1;
    figure(j)
    if j == 1
        R_f_c = linspace(0,0.2,25);
        deg = R_f_c;
    elseif j == 2
        c_tilde_c = linspace(1,0.8,25);
        deg = c_tilde_c;
    elseif j == 3
        c_lyte = linspace(1,0.8,25);
        deg = c_lyte;
    elseif j == 4
        R_f_a = linspace(0,0.2,25);
        deg = R_f_a;
    elseif j == 5
        c_tilde_a = linspace(1,0.8,25);
        deg = c_tilde_a;
    end
   

    %do fro small and large overpotential
    W_charge = solve_W(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, mures_charge_vd_c(1), mures_charge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_discharge = solve_W(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, mures_discharge_vd_c(1), mures_discharge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_charge_exact = solve_W_exact(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, mures_charge_vd_c(1), mures_charge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    W_discharge_exact = solve_W_exact(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, mures_discharge_vd_c(1), mures_discharge_vd_a(1), k0_c, k0_a, lambda_c, lambda_a);
    voltage_charge = solve_W_exact_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_charge_v_c, c_charge_v_a, k0_c, k0_a, lambda_c, lambda_a,1);
    voltage_discharge = solve_W_exact_current(R_f_c, c_tilde_c, c_lyte, R_f_a, c_tilde_a, params_nondim, ...
        c_discharge_v_c, c_discharge_v_a, k0_c, k0_a, lambda_c, lambda_a,-1);
    if j == 5
        plot(deg, W_charge_exact, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', '--')
        hold on
        box on
        plot(deg, W_discharge_exact, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', '--')
        plot(deg, voltage_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', ':')
        plot(deg, voltage_discharge, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', ':')
        legend('Voltage Discharge', 'Voltage Charge', 'Current Discharge', 'Current Charge', 'Fontsize', 16, 'Location', 'southwest')
        legend boxoff
        set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')

    % elseif j == 1
    %     plot(deg, W_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', '--')
    %     hold on
    %     box on
    %     plot(deg, W_discharge, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', '--')
    %     plot(deg, voltage_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', ':', 'HandleVisibility', 'off')
    %     plot(deg, voltage_discharge, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', ':', 'HandleVisibility', 'off')
    %     legend('Charge', 'Discharge', 'Fontsize', 16, 'Location', 'southwest')
    %     legend boxoff
    %     set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')

    else 
        plot(deg, W_charge_exact, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', '--', 'HandleVisibility', 'off')
        hold on
        box on
        plot(deg, W_discharge_exact, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', '--', 'HandleVisibility', 'off')
        plot(deg, voltage_charge, 'LineWidth',2,'Color', map_blue_less(end-3,:), 'LineStyle', ':', 'HandleVisibility', 'off')
        plot(deg, voltage_discharge, 'LineWidth',2,'Color', map_blue_less(end,:), 'LineStyle', ':', 'HandleVisibility', 'off')
        set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')

    end
% 
%     if j == 3
%         ylim([0.5, 1.5])
%     end

    set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
    export_fig(sprintf('fig_sensitivity%d.png', j), '-transparent', '-m3', '-r300')

        
end


end



function out = nearZero(x, tol)
out = (abs(x) < tol);
end


