clear;clc;
params = {};

maxNumCompThreads(1);

%nondimesionalization parameters;
params.N_A = 6.02e23; %/mol
params.kBT = 4.11e-21; %J
params.F = 96485; %C/mol

params.R = 0.1e-6; %m, particle size
params.lambda = 3.78*params.kBT; %units of kBT. Zhang et al 2021

params.cmax = 1.3793e28; % maximum concentration in units of charge; sites/m^3
%params.c0 = 0.24; %nondimensionalized by cmax
%params.max_cap = 0.6; %max capacity fraction cycled (0.47/0.42)
params.c0 = 0.3; %nondimensionalized by cmax
params.max_cap = 0.8; %max capacity fraction cycled (0.47/0.42)
params.D_c = 1e-18; %solid diffusivity, m^2/s of standard phase
%params.D_c = 1e-12; %solid diffusivity, m^2/s of standard phase

params.num_cycles = 6;
params.k = 8; % reaction rate (A/m^2); CIET parameters from Zhang et al.
params.k = params.k*6.241e18; %units of charge/(s*m^2)

params.m = 200; %number of discretizations in secondary particle
%discretization
params.dx = params.R/params.m; %secondary particle

params.control = 1; %C-rate = 0.1. -0.1 is charge. +0.1 is discharge
%since we're only really filling 56%, this is closer to 0.1C

disc_radius = linspace(1,0,params.m+1);
params.vf = disc_radius(1:end-1).^3-disc_radius(2:end).^3;

params_nondim = {};
%constants
params_nondim.d = 2e-10;
params_nondim.e = 1.602e-19;
params_nondim.kBT = 298*1.381e-23;
params_nondim.epsilon0 = 8.854e-12;

params_nondim.v = 1;
%params_nondim.v = 0.1;
params_nondim.R_f = 0;
params_nondim.c_lyte = 1;

params_nondim.x = 0.5; %percentage of nickel
params_nondim.y = 0.3; %percentage of maganese
params_nondim.z = 0.2; %percentage of cobalt
params_nondim.e = 1.602e-19;
params_nondim.kBT = 298*1.381e-23;
params_nondim.muR_ref = 0;
params_nondim.muR_ref = -mu_c(params.c0, params_nondim);
params_nondim.L_ref = params.R;
params_nondim.R = params.R/params_nondim.L_ref; %we set this to reference size
params_nondim.D_ref = params.D_c;
params_nondim.D_c = params.D_c/params_nondim.D_ref;
params_nondim.tau = params_nondim.L_ref^2/params_nondim.D_ref; % s, time constant
%decrease time constant to increase the ode solution resolution
params_nondim.dx = params.dx/params_nondim.L_ref;
params_nondim.k = params.k*params_nondim.tau/(params_nondim.L_ref*params.cmax);
params_nondim.control_dim = 4/3*pi*params.R^3*params.cmax*params.control/3600; %units of sites/s, dimensionalized
params_nondim.control = params_nondim.control_dim.*params_nondim.tau/(params.cmax*params_nondim.L_ref^3); %nondimensionalize
params_nondim.control_abs = abs(params_nondim.control);
params_nondim.lambda = params.lambda/params.kBT;

params.c_s_ind = params.m;
params.mu_res_ind = params.m+1;

%voltage cutoffs for graphite
params_nondim.mu_end_low = 0; %0 is 0.26, 9 is 0.36
params_nondim.mu_end_up = 31; %0.89 voltage--31

sol_total1 = HPPC_function(params, params_nondim);
params_nondim.v = 0.9;
% sol_total2 = HPPC_function(params, params_nondim);
params_nondim.v = 1;
params_nondim.R_f = 0.02;
sol_total3 = HPPC_function(params, params_nondim);
params_nondim.R_f = 0;
params_nondim.c_lyte = 0.9;
% sol_total4 = HPPC_function(params, params_nondim);

map_purple_less = [224,236,244
191,211,230
158,188,218
140,150,198
140,107,177
136,65,157
129,15,124
77,0,75]/255;

%create a function that goes between the values of 

params.num_cycles = 2;
params_nondim.v = 1;
params_nondim.R_f = 0;
params_nondim.c_lyte = 1;
sol_total1_volt = HPPC_function_voltage(params, params_nondim);
params_nondim.v = 0.9;
params_nondim.R_f = 0.02;
params_nondim.c_lyte = 0.9;
sol_total3_volt = HPPC_function_voltage(params, params_nondim);


%%% with current stuff

figure(1)
colororder({'k','k'})
yyaxis left
xsmall = 0.07; xlarge = 0.077;
idx = find(sol_total1.x < xsmall | sol_total1.x > xlarge);
plot(res_xvals(sol_total1.x(idx), xsmall, xlarge), -params_nondim.kBT*(sol_total1.y(end,idx)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'Color', map_purple_less(end,:));
hold on
idx = find(sol_total3.x < xsmall | sol_total3.x > xlarge);
plot(res_xvals(sol_total3.x(idx), xsmall, xlarge), -params_nondim.kBT*(sol_total3.y(end,idx)-params_nondim.muR_ref)/params_nondim.e, 'LineStyle', '-', 'LineWidth', 2, 'Color', map_purple_less(end-4,:))
xlim([0.0685, 0.0715])
ylim([4.237, 4.256])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
yyaxis right
xx = linspace(0.0682, 0.0785,500);
idx = find(xx < xsmall | xx > xlarge);
plot(res_xvals(xx(idx), xsmall, xlarge), return_current(xx(idx), params.control), 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')

ylim([-1.2,1.2])
xlim([0.0685, 0.0715])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')


legend({['Voltage' char(10) '(No Degradation)'], ['Voltage' char(10) '(With Degradation)'], 'Current'}, 'FontSize', 16, 'Location', 'south')
legend boxoff
export_fig fig1.png -transparent -m3 -r300




figure(2)
colororder({'k','k'})
yyaxis left
plot(sol_total1.x, -params_nondim.kBT*(sol_total1.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
hold on
%plot(sol_total2.x, -params_nondim.kBT*(sol_total2.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
plot(sol_total3.x, -params_nondim.kBT*(sol_total3.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineStyle', '-', 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
%plot(sol_total4.x, -params_nondim.kBT*(sol_total4.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
xlim([0.06902, 0.0692])
ylim([4.237, 4.247])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% yyaxis right
% xx = linspace(0.06902, 0.0692);
% plot(xx, return_current(xx, params.control), 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')
% ylim([-0.2,1.2])
% xlim([0.06902, 0.0692])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')
export_fig fig2.png -transparent -m3 -r300


figure(3)
colororder({'k','k'})
yyaxis left
plot(sol_total1.x, -params_nondim.kBT*(sol_total1.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
hold on
%plot(sol_total2.x, -params_nondim.kBT*(sol_total2.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
plot(sol_total3.x, -params_nondim.kBT*(sol_total3.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineStyle', '-', 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
%plot(sol_total4.x, -params_nondim.kBT*(sol_total4.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
xlim([0.07775, 0.07795])
ylim([4.249, 4.256])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')
% yyaxis right
% xx = linspace(0.07775, 0.07795);
% plot(xx, return_current(xx, params.control), 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')
% ylim([-1.2,0.2])
% xlim([0.07775, 0.07795])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')

export_fig fig3.png -transparent -m3 -r300




figure(4)

colororder({'k','k'})
yyaxis left
xsmall = 0.07; xlarge = 0.077;
idx1 = find(sol_total1_volt.x < xsmall | sol_total1_volt.x > xlarge);
plot(res_xvals(sol_total1_volt.x(idx1), xsmall, xlarge), sol_total1_volt.y(end,idx1), 'LineWidth', 2, 'Color', map_purple_less(end,:))
hold on
%plot(sol_total2.x, -params_nondim.kBT*(sol_total2.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'Color', map_purple_less(end,:))
idx = find(sol_total3_volt.x < xsmall | sol_total3_volt.x > xlarge);
plot(res_xvals(sol_total3_volt.x(idx), xsmall, xlarge), sol_total3_volt.y(end,idx), 'LineStyle', '-', 'LineWidth', 2, 'Color', map_purple_less(end-4,:))
%plot(sol_total4.x, -params_nondim.kBT*(sol_total4.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'Color', map_purple_less(end-4,:))
xlim([0.0686, 0.0715])
ylim([-23, 23])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
yyaxis right
ylim([4.12, 4.36])
plot(res_xvals(sol_total1_volt.x(idx1), xsmall, xlarge), -params_nondim.kBT*(sol_total1_volt.y(end-1,idx1)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')
xlim([0.0686, 0.0715])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')


legend({['Current' char(10) '(No Degradation)'], ['Current' char(10) '(With Degradation)'], 'Voltage'}, 'FontSize', 16, 'Location', 'south')
legend boxoff
export_fig fig1_volt.png -transparent -m3 -r300


figure(5)
colororder({'k','k'})
yyaxis left
plot(sol_total1_volt.x, sol_total1_volt.y(end,:), 'LineWidth', 4, 'Color', map_purple_less(end,:))
hold on
%plot(sol_total2.x, -params_nondim.kBT*(sol_total2.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
plot(sol_total3_volt.x, sol_total3_volt.y(end,:), 'LineStyle', '-', 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
%plot(sol_total4.x, -params_nondim.kBT*(sol_total4.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
xlim([0.06902, 0.0692])
ylim([-23, -7])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% yyaxis right
% plot(sol_total1_volt.x, -params_nondim.kBT*(sol_total1_volt.y(end-1,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')
% ylim([4.20, 4.38])
% xlim([0.06902, 0.0692])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')
export_fig fig2_volt.png -transparent -m3 -r300


figure(6)
colororder({'k','k'})
yyaxis left
plot(sol_total1_volt.x, sol_total1_volt.y(end,:), 'LineWidth', 4, 'Color', map_purple_less(end,:))
hold on
%plot(sol_total2.x, -params_nondim.kBT*(sol_total2.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end,:))
plot(sol_total3_volt.x, sol_total3_volt.y(end,:), 'LineStyle', '-', 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
%plot(sol_total4.x, -params_nondim.kBT*(sol_total4.y(end,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 4, 'Color', map_purple_less(end-4,:))
xlim([0.07775, 0.07795])
ylim([7, 23])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')
% yyaxis right
% xx = linspace(0.07775, 0.07795);
% plot(sol_total1_volt.x, -params_nondim.kBT*(sol_total1_volt.y(end-1,:)-params_nondim.muR_ref)/params_nondim.e, 'LineWidth', 2, 'LineStyle',':', 'Color', 'k')
% xlim([0.07775, 0.07795])
% ylim([4.12, 4.28])
set(gca,'fontsize',12,'LineWidth',4,'FontName','Helvetica')

export_fig fig3_volt.png -transparent -m3 -r300





% xlabel('Time (h)')
% ylabel('Voltage (V)')
%params_nondim.v


function sol_total = HPPC_function(params, params_nondim)

% set up volume rate 
M = spdiags(ones(params.c_s_ind,1),0,params.m, params.m);
M(params.mu_res_ind, params.mu_res_ind) = 0;
M(params.mu_res_ind, 1:params.c_s_ind) = params.vf;
% also add the v_params for current control
%M(params.mu_res_ind, params.v_ind+1:params.delta_ind) = params.vf;

%find initial condition
c_s0 = [ones(params.m,1)*params.c0];

%options = optimset('TolFun', 1e-10, 'TolX', 1e-10, 'Display', 'none');
options = optimset('TolFun', 1e-8, 'TolX', 1e-8);

%how many cycles do we want to do?
sol_total = {};
sol_total.x = [];
sol_total.y = [];
t_charge_end = [];
t_discharge_end = [];
t_end = 0;


for i = 1:params.num_cycles
    
    i
    %step 1: charge
    
    params_nondim.control = -params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_up;
    params_nondim.cap_max = params.c0+params.max_cap;
    %'RelTol', 1e-7, 'AbsTol', 1e-7, 
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), -1, [], [], [], [], [], [], [], options);
 %   -params_nondim.kBT*(sol_total3.y(end,:)-params_nondim.muR_ref)/params_nondim.e,
   % params_nondim.muR_ref-params_nondim.e/params_nondim.kBT*
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(360*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    
    %write down old variables
    t_charge_end = [t_charge_end, sol.x(end)];
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
    
    params_nondim.control = 0;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    %opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(6*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    t_discharge_end = [t_discharge_end, sol.x(end)];
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

    params_nondim.control = params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_up;
    params_nondim.cap_max = params.c0+params.max_cap;
    %'RelTol', 1e-7, 'AbsTol', 1e-7, 
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), -1, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(360*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    
    %write down old variables
    t_charge_end = [t_charge_end, sol.x(end)];
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

    params_nondim.control = 0;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    %opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(6*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    t_discharge_end = [t_discharge_end, sol.x(end)];
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
   
    params_nondim.control = -0.1*params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    %opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    t_discharge_end = [t_discharge_end, sol.x(end)];
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

end

sol.x = sol.x*params_nondim.tau;

end


function sol_total = HPPC_function_voltage(params, params_nondim);

% set up volume rate 
M = spdiags(ones(params.c_s_ind,1),0,params.m, params.m);
M(params.mu_res_ind, params.mu_res_ind) = 0;
M(params.mu_res_ind, 1:params.c_s_ind) = params.vf;
% also add the v_params for current control
%M(params.mu_res_ind, params.v_ind+1:params.delta_ind) = params.vf;

%find initial condition
c_s0 = [ones(params.m,1)*params.c0];

%options = optimset('TolFun', 1e-10, 'TolX', 1e-10, 'Display', 'none');
options = optimset('TolFun', 1e-8, 'TolX', 1e-8);

%how many cycles do we want to do?
sol_total = {};
sol_total.x = [];
sol_total.y = [];
t_charge_end = [];
t_discharge_end = [];
t_end = 0;


for i = 1:params.num_cycles
    
    i
    %step 1: charge
    pulse = -0.01*(params_nondim.e/params_nondim.kBT);
    
    params_nondim.control = -params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_up;
    params_nondim.cap_max = params.c0+params.max_cap;
    %'RelTol', 1e-7, 'AbsTol', 1e-7, 
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), -1, [], [], [], [], [], [], [], options);
   % IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
    %IC = [c_s0; real(IC_0)];
    IC = c_s0;
    sol = ode15s(@(t,f) dfdt_voltage(t, f, params, params_nondim, mu_c(params.vf*c_s0, params_nondim)+pulse), [0, 0.6/(360*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    dy = dfdt_voltage_plot_only(x, y, params, params_nondim, mu_c(params.vf*c_s0, params_nondim)+pulse);
    curr = params.vf*dy;
    
    %write down old variables
    t_charge_end = [t_charge_end, sol.x(end)];
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, [y; (mu_c(params.vf*c_s0, params_nondim)+pulse*10)*ones(1,80); curr]];
   

    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
    
    params_nondim.control = 0;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c(params.vf*c_s0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(6*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    t_discharge_end = [t_discharge_end, sol.x(end)];
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, [y; zeros(1,80)]];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

    params_nondim.control = params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_up;
    params_nondim.cap_max = params.c0+params.max_cap;
    %'RelTol', 1e-7, 'AbsTol', 1e-7, 
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
   % opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), -1, [], [], [], [], [], [], [], options);
   % IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c((i-1)*0.1+params.c0, params_nondim), options);
   % IC = [c_s0; real(IC_0)];
    IC = c_s0;
    sol = ode15s(@(t,f) dfdt_voltage(t, f, params, params_nondim, mu_c(params.vf*c_s0, params_nondim)-pulse), [0, 0.6/(360*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    dy = dfdt_voltage_plot_only(x, y, params, params_nondim, mu_c(params.vf*c_s0, params_nondim)-pulse);
    curr = params.vf*dy;

    
    %write down old variables
    t_charge_end = [t_charge_end, sol.x(end)];
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, [y; (mu_c(params.vf*c_s0, params_nondim)-pulse*10)*ones(1,80); curr]];
    
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

    params_nondim.control = 0;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    %opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c(params.vf*c_s0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(6*params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);

    t_discharge_end = [t_discharge_end, sol.x(end)];

    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, [y; zeros(1,80)]];

    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
   
    params_nondim.control = -0.1*params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;

    %'RelTol', 1e-7, 'AbsTol', 1e-7,
    %opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    opts = odeset('Mass', M, 'RelTol', 1e-7, 'AbsTol', 1e-7, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    %IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 30, [], [], [], [], [], [], [], options);
    IC_0 = fsolve(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), mu_c(params.vf*c_s0, params_nondim), options);
    IC = [c_s0; real(IC_0)];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, 0.6/(params_nondim.control_abs)], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    t_discharge_end = [t_discharge_end, sol.x(end)];
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, [y;  zeros(1,80)]];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);

end

sol.x = sol.x*params_nondim.tau;

end




function [value,isterminal,direction] = bounceEvents(t,y,params,params_nondim)
% %end simulation when certain mures is reached
% %y(1)./(1-y(params.c_s_ind+1))-params_nondim.cap_max;
% value = [params.vf*y(1:params.c_s_ind)-params_nondim.cap_max; 
%     params.vf*y(1:params.c_s_ind)-params.c0];
% % value = [y(params.mu_res_ind)-params_nondim.mu_end;
% %     params.vf*y(1:params.c_s_ind)-params_nondim.cap_max; 
% %     params.vf*y(1:params.c_s_ind)-params.c0];
% % first one detects cutoff voltage, second one is the filling at the 
% isterminal = [1; 1];   % Stop the integration
% direction = [0; 0];   % Negative direction only
% % isterminal = [1; 1; 1];   % Stop the integration
% % direction = [0; 0; 0];   % Negative direction only
value = [params.vf*y(1:params.c_s_ind)-0.95; params.vf*y(1:params.c_s_ind)-0.2];
isterminal = [1; 1];   % Stop the integration
direction = [0; 0];   % Negative direction only
end


function solve_y = initial_condition(IC0, M, cs_0, params, params_nondim)
%solves for initial condition of DAE
y0 = [cs_0; IC0];
y = dfdt(0, y0, params, params_nondim);
solve_y = M(params.mu_res_ind,1:params.c_s_ind)*y(1:params.c_s_ind)-y(params.mu_res_ind);
end


function out = D(c, params_nondim)
% reutrns diffusivity of the "blocked sections using c_tilde"
%out = params_nondim.D_c.*ctilde.*(1-ctilde);
% for concentration, 
%out = params_nondim.D_c;
out = 20000.*10.^(-2319.*c.^10 + 6642.*c.^9 - 5269.*c.^8 - 3319.*c.^7 ...
    + 10038.*c.^6 - 9806.*c.^5 + 5817.*c.^4 - 2286.*c.^3 + 575.3*c.^2 ...
    - 83.16*c-9.292)/1e2/params_nondim.D_ref;
end


function out = mu_c(c, params_nondim)
%here we have the scaled chemical potential summed from the two different
%phases
%OCV from Colclasure 2020 NMC532
OCV = Colclasure_OCV(c);
%the OCV is in units of V, convert to kBT
out = -params_nondim.e/params_nondim.kBT*OCV + params_nondim.muR_ref;
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



function out = R_value(rxn, c, v, R_f, c_lyte, muh, mures, params_nondim)
% c_eff = c./(1-v); %scale by # of available sites
% %c_eff = c;
% eta = (muh - mures);
% alpha = 0.5;
% % rxn = params_nondim.k.*(165*c_eff.^5 - 752.4*c_eff.^4 + 1241*c_eff.^3 ...
% %       - 941.7*c_eff.^2 + 325*c_eff - 35.85).*(exp(-alpha*eta)-exp((1-alpha)*eta));
% rxn = params_nondim.k.*sqrt(c_eff.*(1-c_eff)).*(exp(-alpha*eta)-exp((1-alpha)*eta));
eta = muh-mures;
eta_f = muh - mures - log(c./c_lyte) + rxn*R_f;
i_red = helper_fun(-eta_f, params_nondim.lambda);
i_ox = helper_fun(eta_f, params_nondim.lambda);
out = -params_nondim.k/sqrt(4*pi*params_nondim.lambda)*(v-c).*(c_lyte.*i_red - c.*i_ox) - rxn;
%rxn = params_nondim.k/sqrt(4*pi*params_nondim.lambda)*(1-c).*(i_red - c.*i_ox);
end


function out = R(c, v, R_f, c_lyte, muh, mures, params_nondim)
out = fzero(@(rxn) R_value(rxn, c, v, R_f, c_lyte, muh, mures, params_nondim), 0);
end


function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end


function output = return_current(time, control)
%get current value
idx_pos = (time > 0.0690414) & (time < 0.0691847);
idx_neg = (time > 0.077779) & (time < 0.0779223);
output = time * 0;
output(idx_pos) = control;
output(idx_neg) = -control;
end


function dy = dfdt_voltage(t, y, params, params_nondim, mu_res)
%y = c(1..m); c_gb(1..m); mures
%assuming homogenuous particles and surface reaction
%get flux terms for the particles

%flux condition between grain boundary and particle
%we assume k = k*A area
%lattice model for diffusion in primary particle
radii = linspace(params_nondim.R,0,params.m+1).'; %nondimensional radius
radii_mid = (radii(1:end-1)+radii(2:end))/2;
c_s = y;
muc = mu_c(c_s, params_nondim);

Rxn = R(c_s(1), params_nondim.v, params_nondim.R_f, params_nondim.c_lyte, muc(1), mu_res, params_nondim);

F_D = -diff(D(c_s, params_nondim).*c_s.*muc)/params_nondim.dx;
%get the boundary fluxes
F_bc = -Rxn;   
%we are assuming surface reaction for the big particle
F = [F_bc; F_D; 0];
dy(1:params.c_s_ind,1) = -1./radii_mid.^2.*diff(radii.^2.*F)/params_nondim.dx;
%first term is the same as prevoius
end



function dy = dfdt_voltage_plot_only(t, y, params, params_nondim, mu_res)
% for plotting purpsoes only
%assuming homogenuous particles and surface reaction
%get flux terms for the particles

%flux condition between grain boundary and particle
%we assume k = k*A area
%lattice model for diffusion in primary particle
radii = linspace(params_nondim.R,0,params.m+1).'; %nondimensional radius
radii_mid = (radii(1:end-1)+radii(2:end))/2;
c_s = y;
muc = mu_c(c_s, params_nondim);

Rxn = zeros(size(y,2),1);
for j = 1:size(y,2)
    Rxn(j) = R(c_s(1,j), params_nondim.v, params_nondim.R_f, params_nondim.c_lyte, muc(1,j), mu_res, params_nondim);
end

F_D = -diff(D(c_s, params_nondim).*c_s.*muc)/params_nondim.dx;
%get the boundary fluxes
F_bc = -Rxn;   
%we are assuming surface reaction for the big particle
F = [F_bc.'; F_D; zeros(1,80)];
dy = -1./radii_mid.^2.*diff(radii.^2.*F)/params_nondim.dx;
%first term is the same as prevoius
end



function dy = dfdt(t, y, params, params_nondim)
%first get voltage solution
mu_res = y(params.mu_res_ind);
dy = dfdt_voltage(t, y(1:params.c_s_ind), params, params_nondim, mu_res);

%last term is an algebraic equation
dy(params.mu_res_ind,1) = -params_nondim.control;
end



function h=BreakXAxis(x,y,start,stop,width,lwidth,color)
% Julie Haas, after Michael Robbins
% - assumes dx is same for all data
% - 'width' must be a whole number
% - to change axis, use {a=axis} to get axis, and reset in those units
%test data:  
% x=[1:.5:40]; y=rand(size(x)); start=10;stop=20;width=6; 
% x=.01:.01:10;y=sin(6*x);start=2;stop=3;width=1;
% erase unused data
y(x>start & x<stop)=[];
x(x>start & x<stop)=[];
% map to new xaxis, leaving a space 'width' wide
x2=x;
x2(x2>=stop)=x2(x2>=stop)-(stop-start-width);
h=plot(x2,y,':','Linewidth',lwidth,'Color',color);
ytick=get(gca,'YTick');
t1=text(start+width/2,ytick(1),'//','fontsize',15);
t2=text(start+width/2,ytick(max(length(ytick))),'//','fontsize',15);
% For y-axis breaks, use set(t1,'rotation',270);
% remap tick marks, and 'erase' them in the gap
xtick=get(gca,'XTick');
dtick=xtick(2)-xtick(1);
gap=floor(width/dtick);
last=max(xtick(xtick<=start));          % last tick mark in LH dataset
next=min(xtick(xtick>=(last+dtick*(1+gap))));   % first tick mark within RH dataset
offset=size(x2(x2>last&x2<next),2)*(x(2)-x(1));
for i=1:sum(xtick>(last+gap))
    xtick(find(xtick==last)+i+gap)=stop+offset+dtick*(i-1);
end
    
for i=1:length(xtick)
    if xtick(i)>last&xtick(i)<next
        xticklabel{i}=sprintf('%d',[]);
    else
        xticklabel{i}=sprintf('%d',xtick(i));
    end
end;
set(gca,'xticklabel',xticklabel);
end


 function xvals = res_xvals(xvals, xsmall, xlarge)
 idx = find(xvals > xlarge);
 xvals(idx) = xvals(idx) - xlarge + xsmall; 
 end

%%% plotting useful lines
