params = {};

%%%create script for running constant current cycling experiments with the defect model

%nondimesionalization parameters;
params.N_A = 6.02e23; %/mol
params.kBT = 4.11e-21; %J
params.F = 96485; %C/mol

params.R = 0.1e-6; %m, particle size
params.lambda = 3.78*params.kBT; %units of kBT. Zhang et al 2022

params.cmax = 1.3793e28; % maximum concentration in units of charge; sites/m^3
params.c0 = 0.4; %nondimensionalized by cmax
%params.max_cap = 0.8; %max capacity fraction cycled (0.47/0.42)
params.max_cap = 2; %max capacity fraction cycled (0.47/0.42)
params.v0 = 0.02; %max vacancy concentration
params.D_c = 1e-16; %solid diffusivity, m^2/s of standard phase
params.D_v = 5e-24; %diffusivity of nickel
params.phi0_oxy = 4.4; %voltage cutoff for formation of oxygen va
%cutoff for 4.8 V: 0.2601
%cutoff for 4.6 V: 0.2620
%cutoff for 4.4 V: 0.2838
%cutoff for 4.2 V: 0.3629

params.num_cycles = 10;
params.k = 8; % reaction rate (A/m^2); CIET parameters from Zhang et al.
params.k = params.k*6.241e18; %units of charge/(s*m^2)
params.k_oxy = 1e-5; % reaction rate (A/m^2)
params.k_oxy = params.k_oxy*6.241e18; %units of charge/(s*m^2)

params.m = 200; %number of discretizations in secondary particle
%discretization
params.dx = params.R/params.m; %secondary particle

params.control = 1; %C-rate = 0.1. -0.1 is charge. +0.1 is discharge
params.control = params.control * 0.53; %since we assume 53% is the max capacity

params.mu_res_low = 4.5;
params.mu_res_high = 3.75;
params.cap_res_high = 0.9;

disc_radius = linspace(1,0,params.m+1);
params.vf = disc_radius(1:end-1).^3-disc_radius(2:end).^3;

params_nondim = {};
params_nondim.x = 0.5; %percentage of nickel
params_nondim.y = 0.3; %percentage of maganese
params_nondim.z = 0.2; %percentage of cobalt
params_nondim.e = 1.602e-19;
params_nondim.kBT = 298*1.381e-23;
params_nondim.muR_ref = 0;
params_nondim.muR_ref = -mu_c(params.c0, 0, params_nondim);
params_nondim.L_ref = params.R;
params_nondim.R = params.R/params_nondim.L_ref; %we set this to reference size
params_nondim.D_ref = params.D_c;
params_nondim.D_c = params.D_c/params_nondim.D_ref;
params_nondim.D_v = params.D_v/params_nondim.D_ref;
params_nondim.tau = params_nondim.L_ref^2/params_nondim.D_ref; % s, time constant
%decrease time constant to increase the ode solution resolution
params_nondim.dx = params.dx/params_nondim.L_ref;
params_nondim.k = params.k*params_nondim.tau/(params_nondim.L_ref*params.cmax);
params_nondim.k_oxy = params.k_oxy*params_nondim.tau/(params_nondim.L_ref*params.cmax);
params_nondim.control_dim = 4/3*pi*params.R^3*params.cmax*params.control/3600; %units of sites/s, dimensionalized
params_nondim.control = params_nondim.control_dim.*params_nondim.tau/(params.cmax*params_nondim.L_ref^3); %nondimensionalize
params_nondim.phi0_oxy = -params.phi0_oxy*params_nondim.e/params_nondim.kBT + params_nondim.muR_ref;
params_nondim.control_abs = abs(params_nondim.control);
params_nondim.lambda = params.lambda/params.kBT;

params.c_s_ind = params.m;
params.v_ind = 2*params.m;
params.oxy_ind = 2*params.m+1;
params.mu_res_ind = 2*params.m+2;

c_v = load('Ni0.5Mn0.3Co0.2_c_v.txt');
params_nondim.c_av = c_v(:,1);
params_nondim.v_av = c_v(:,2);
%these only contain enthalpic interactions, entropic are calculated
%analytically to avoid 0s
params_nondim.muc = load('Ni0.5Mn0.3Co0.2_muc_ent.txt');
params_nondim.muv = load('Ni0.5Mn0.3Co0.2_muv_ent.txt');

%voltage cutoffs for graphite
params_nondim.mu_end_low = -params.mu_res_low*params_nondim.e/params_nondim.kBT + params_nondim.muR_ref; %0 is 0.26, 9 is 0.36
params_nondim.mu_end_up = -params.mu_res_high*params_nondim.e/params_nondim.kBT + params_nondim.muR_ref; %0.89 voltage--31

% set up volume rate 
M = diag(ones(params.oxy_ind,1));
M(params.mu_res_ind, 1:params.c_s_ind) = params.vf;
M(params.mu_res_ind, params.mu_res_ind) = 0;

%create interpolating function
params_nondim.muv_interp = griddedInterpolant({params_nondim.c_av, params_nondim.v_av}, params_nondim.muv.', 'linear', 'spline');

%find initial condition
c_s0 = [ones(params.m,1)*params.c0];
v0 = [ones(params.m,1)*params.v0; 1]; %adding the oxygen concentratoin in the v array

options = optimset('TolFun', 1e-10, 'TolX', 1e-10, 'Display', 'none');
opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8);

%how many cycles do we want to do?
sol_total = {};
sol_total.x = [];
sol_total.y = [];
t_end = 0;

%also store max, min capacities
max_capacities = [];
min_capacities = [];

for i = 1:params.num_cycles
    
    i
    
    params_nondim.control = -params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_up;
    params_nondim.cap_max = params.c0+params.max_cap;
    opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), -1, [], [], [], [], [], [], [], options);
    IC = [c_s0; v0; IC_0];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, params.max_cap/params_nondim.control_abs], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
    v0 = sol_total.y(params.c_s_ind+1:params.oxy_ind,end);
    max_capacities = [max_capacities, params.vf*sol_total.y(1:params.c_s_ind,end)/params.m];
    
    params_nondim.control = params_nondim.control_abs;
    params_nondim.mu_end = params_nondim.mu_end_low;
    params_nondim.cap_max = params.c0;
    opts = odeset('Mass', M, 'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params, params_nondim));
    IC_0 = fmincon(@(inputs) initial_condition(inputs, M, c_s0, params, params_nondim), 15, [], [], [], [], [], [], [], options);
    IC = [c_s0; v0; IC_0];
    sol = ode15s(@(t,f) dfdt(t, f, params, params_nondim), [0, params.max_cap/params_nondim.control_abs], IC, opts);
    x = linspace(0, sol.x(end), 80);
    y = deval(sol, x);
    
    %write down old variables
    sol_total.x = [sol_total.x, t_end + x];
    sol_total.y = [sol_total.y, y];
    
    %write down "end" variables
    t_end = sol_total.x(end);
    c_s0 = sol_total.y(1:params.c_s_ind,end);
    v0 = sol_total.y(params.c_s_ind+1:params.oxy_ind,end);
    min_capacities = [min_capacities, params.vf*sol_total.y(1:params.c_s_ind,end)/params.m];


end

sol.x = sol.x*params_nondim.tau;

save('Ni' + string(params_nondim.x) + '_' + string(params_nondim.control_abs) + 'C_5_cycle_high.mat')

% h1 = animatedline();
% axis([0,1,0,1]);
% spacings = linspace(0,1,params.m);
% 
% xlabel('Length (nondim)')
% ylabel('Solid Concentration')
% 
% for k = ceil(linspace(0,sol_total.x(end),400))
%     clearpoints(h1);
%     if k == 0
%         k = 1;
%     end
%     addpoints(h1,spacings,sol_total.y(1:params.c_s_ind,k));
%     drawnow
% end
% 
% figure(2)
% h2 = animatedline();
% axis([0,1,0,1]);
% 
% xlabel('Length (nondim)')
% ylabel('Disorder Concentration')
% 
% for k = ceil(linspace(0,sol_total.x(end),400))
%     clearpoints(h2);
%     if k == 0
%         k = 1;
%     end
%     addpoints(h2,spacings,1-sol_total.y(params.c_s_ind+1:params.v_ind,k));
%     drawnow
% end
% 
% save('Ni' + string(params_nondim.x) + '_' + string(params_nondim.control_abs)...
%     + 'C.mat')


%overpotential for disorder reaction;
eta_v = -(sol_total.y(params.mu_res_ind,:) - params_nondim.phi0_oxy);

zci = @(v) find(v(:).*circshift(v(:),[-1 0]) <= 0);                    % Returns Approximate Zero-Crossing Indices Of Argument Vector
dy = zci(eta_v);                                                            % Indices of Approximate Zero-Crossings
x0 = zeros(size(dy));
for k1 = 1:size(dy,1)
    b = [[1;1] [sol_total.x(dy(k1)); sol_total.x(dy(k1)+1)]]\[eta_v(dy(k1)); eta_v(dy(k1)+1)];      % Linear Fit Near Zero-Crossings
    x0(k1) = -b(1)/b(2);                                                % Interpolate ‘Exact’ Zero Crossing
    mb(:,k1) = b;                                                       % Store Parameter Estimates (Optional)
end


figure(3)
%plot amount of disorder
dis = sol_total.y(params.c_s_ind+1,:);
plot(sol_total.x, dis)
hold on
xlabel('Time (s)')
ylabel('Disordered Concentration')
ylim([min(dis)-0.00001, max(dis)+0.00001])
for i = 1:length(x0)
    plot(x0(i)*ones(200,1), linspace(0,1,200), '--r')
end


figure(4)
%plot overpotential of eta
plot(sol_total.x, eta_v)
hold on
xlabel('Time (s)')
ylabel('\eta_v')
ylim([min(eta_v)-4,max(eta_v)+4])
for i = 1:length(x0)
    plot(x0(i)*ones(200,1), linspace(-100,100,200), '--r')
end
plot(sol_total.x, 0*sol_total.x, '--k')



function eps = epsilon(c, v, delta, params_nondim)
%create function for nonlinear epsilon
alpha = 0.03*c + 0.266*(params_nondim.x-v) + 0.544*params_nondim.y + 0.508*params_nondim.z + 2*2.75*delta;
%molecular mass
MM = 6.941*c + 58.693*(params_nondim.x-v) + 54.938*params_nondim.y + 58.933*params_nondim.z + delta*2*15.999;
%number density
%density of NMC~ 2e6 per A^3
prefac = 6.02e23*2.3e6/(1e30);
numdensity = prefac./MM;
eps = (-3-8*alpha.*numdensity*pi)./(-3+4*alpha.*numdensity*pi);
end



function [value,isterminal,direction] = bounceEvents(t,y,params,params_nondim)
%end simulation when certain mures is reached
%y(1)./(1-y(params.c_s_ind+1))-params_nondim.cap_max;
value = [y(params.mu_res_ind)-params_nondim.mu_end_low;
    y(params.mu_res_ind)-params_nondim.mu_end_up;
    params.vf*y(1:params.c_s_ind) - params.cap_res_high];
% first one detects cutoff voltage, second one is the filling at the 
isterminal = [1; 1; 1];   % Stop the integration
direction = [0; 0; 0];   % Negative direction only
end


function solve_y = initial_condition(IC0, M, cs_0, params, params_nondim)
%solves for initial condition of DAE
y0 = [cs_0; params.v0*ones(size(cs_0)); 1; IC0];
y = dfdt(0, y0, params, params_nondim);
solve_y = abs(M(params.mu_res_ind,1:params.c_s_ind)*y(1:params.c_s_ind)-y(params.mu_res_ind));
end


function out = D(c, v, params_nondim)
% reutrns diffusivity of the "blocked sections using c_tilde"
%out = params_nondim.D_c.*ctilde.*(1-ctilde);
% for concentration, 
out = params_nondim.D_c;
% out = 20000.*10.^(-2319.*c.^10 + 6642.*c.^9 - 5269.*c.^8 - 3319.*c.^7 ...
%     + 10038.*c.^6 - 9806.*c.^5 + 5817.*c.^4 - 2286.*c.^3 + 575.3*c.^2 ...
%     - 83.16*c-9.292)/1e4/params_nondim.D_ref;
end


function out = mu_c(c, v, params_nondim)
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


function out = mu_v(c, v, delta, params_nondim)
x = params_nondim.x;
%plut in negative value of entropy because v is negative
%muv_ent = real(log(v./(1-c-v))+log(v./(x-v))+1e6.*heaviside(1-c-v).*(1-c-v)+1e6.*heaviside(x-v).*(x-v));
%muv_ent = log(v./(1-c-v))+log(v./(x-v));% - params_nondim.kappa*sphere_wet(v, r_vec, params_nondim.dx);
muv_ent = log(v./(1-c-v))+log(v./(x-v));% - params_nondim.kappa*sphere_wet(v, r_vec, params_nondim.dx);
%get the enthalpic terms
out = (params_nondim.e/params_nondim.kBT*params_nondim.muv_interp(c, v))./epsilon(c, v, delta, params_nondim).*epsilon(c, v, 1, params_nondim) +200;
%out = (params_nondim.e/params_nondim.kBT*params_nondim.muv_interp(c, v)+200)./delta;
%switch to directino of chemical ptoential and not voltage
%the enthalpic terms are in units of V, convert to kBT
out = (out + muv_ent);
%out = muv_ent;
%out = muv_ent;
end



function rxn = R(c, v, muh, mures, params_nondim)
% c_eff = c./(1-v); %scale by # of available sites
% %c_eff = c;
% eta = (muh - mures);
% alpha = 0.5;
% % rxn = params_nondim.k.*(165*c_eff.^5 - 752.4*c_eff.^4 + 1241*c_eff.^3 ...
% %       - 941.7*c_eff.^2 + 325*c_eff - 35.85).*(exp(-alpha*eta)-exp((1-alpha)*eta));
% rxn = params_nondim.k.*sqrt(c_eff.*(1-c_eff)).*(exp(-alpha*eta)-exp((1-alpha)*eta));
eta = mures-muh;
eta_f = mures + log(c) - muh;
i_red = helper_fun(eta_f, params_nondim.lambda);
i_ox = helper_fun(-eta_f, params_nondim.lambda);
rxn = params_nondim.k/sqrt(4*pi*params_nondim.lambda)*(1-c-v).*(i_red - c.*i_ox);
end


function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end


function rxn = R_v(mures, oxy, params_nondim)
eta = (mures-params_nondim.phi0_oxy);
%this extra negative comes from fucking up a sign somewhere... it's also
%fucked up in the intercalation, but it got carried over to the algebraic
%control equation
rxn = params_nondim.k_oxy.*oxy.*exp(-eta);
end


function dy = dfdt(t, y, params, params_nondim)
%y = c(1..m); c_gb(1..m); mures
%assuming homogenuous particles and surface reaction
%get flux terms for the particles

%flux condition between grain boundary and particle
%we assume k = k*A area
%lattice model for diffusion in primary particle
radii = linspace(params_nondim.R,0,params.m+1).'; %nondimensional radius
radii_mid = (radii(1:end-1)+radii(2:end))/2;

c_s = y(1:params.c_s_ind);
v = y(params.c_s_ind+1:params.v_ind);
oxy = y(params.oxy_ind);
muc = mu_c(c_s, v, params_nondim);
%append boundary condition to muv, c_oxy constant
muv = mu_v(c_s, v, 1, params_nondim);
%we use a "fixed" bc of 50% full outside. this needs to be changed
mu_res = y(params.mu_res_ind);

Rxn = R(c_s(1), v(1), muc(1), mu_res, params_nondim);
%boundary condition oxygen reaction
Rxn_v = R_v(mu_res, oxy, params_nondim);

F_D = -diff(D(c_s, v, params_nondim).*c_s.*(1-v).*muc)/params_nondim.dx;
F_v = -diff(params_nondim.D_v.*v.*muv)/params_nondim.dx;
%get the boundary fluxes
F_bc = -Rxn;  
F_v_bc = -(-Rxn_v);
%we are assuming surface reaction for the big particle
F = [F_bc; F_D; 0];
F_v = [F_v_bc; F_v; 0];
dy(1:params.c_s_ind,1) = -1./radii_mid.^2.*diff(radii.^2.*F)/params_nondim.dx;
dy(params.c_s_ind+1:params.v_ind,1) = -1./radii_mid.^2.*diff(radii.^2.*F_v)/params_nondim.dx;
dy(params.oxy_ind,1) = -Rxn_v*(4/3*pi*(radii(1)^3-radii(2)^3))./params_nondim.dx; 
%first term is the same as prevoius

%last term is an algebraic equation
dy(params.mu_res_ind,1) = -params_nondim.control;
end
