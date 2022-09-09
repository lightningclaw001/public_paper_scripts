clc;clearvars

maxNumCompThreads(1);

%create input params
params = {};
%can be set to current or voltage for CC/CV
params.protocol = "current";
%params.protocol = "voltage";
%control value of current or voltage
params.control_value = 0.4;
%params.control_value = 0;
%set reaction mode, choose between linear and BV right now
%params.rxn_method = "BV";
params.rxn_method = "CIET";
%params.rxn_method = "Marcus";
%set final time
params.final_t = 2;
%set radius dicretization
%initial concentration
params.c0 = 0.45;
params.delta_c = 0.002;
params.delta_R = 10e-9;
params.c_array = params.delta_c+0.35:params.delta_c:1-params.delta_c;
params.R_array = 10e-9:params.delta_R:300e-9;
%nondimensionalize
%params.num_cycles = 24;
%params.num_cycles = 1;
params.num_cycles = 100;
params.R_array = params.R_array/100e-9;
%set 0.1 as the reference concentration
params.muR_ref = 0;
params.e = 1.602e-19;
params.kBT = 298*1.381e-23;
params.lambda = 3.78; %units of kBT. Zhang et al 2021
params.muR_ref = -mu_c(params.c0, params);
%find the reference OCV value for the 
params.V_res = -params.e/params.kBT*4.1 + params.muR_ref;
params.delta_r = params.delta_R/100e-9;
%we need to recenter everything at the middle reaction rate
[c_grid, R_grid] = meshgrid(params.c_array, params.R_array);
%flip array so we can get opposite orderings
params.c_grid = c_grid';
params.R_grid = R_grid';
%set down to 2d array
%first type is ordered (c,T), which is also default ordering
params.c_grid_1 = params.c_grid(:);
params.R_grid_1 = params.R_grid(:);
%reprocess grid to be stacked
params.N = length(params.c_array);
params.M = length(params.R_array);
params.NM = params.N*params.M;
params.sigma_c = 0.02;
%params.sigma_R = 10e-9;
params.sigma_R = 100e-9;
params.avg_R = 100e-9;
%params.beta = 50;
params.beta = 1000;
%params.beta = 0.2;
params.avg_R = params.avg_R/100e-9;
params.sigma_R = params.sigma_R/100e-9;
%params.order = 1;
params.order = 10;

params.cutoff = 2;

%set voltage cutoffs
%if we're using voltage cutoffs
params.voltage_cutoffs = true(1);
params.upper_volt = -params.e/params.kBT*4.1 + params.muR_ref;
params.lower_volt = -params.e/params.kBT*3.7 + params.muR_ref;
params.cap_res_high = 0.95;
params.cap_res_low = 0.2;

%set reaction parameters
rxn_params = {};
%interaction paraemter for LFP
%diffusive coefficeint D0 kBT/N_t
rxn_params.D0 = 5e-2;
%rxn_params.D0 = 1e-2;
%rxn_params.k0 = 25/50;
rxn_params.k0 = 10;
% rxn_params.k0_res = 0.1e-4
rxn_params.k0_res = 0.002;
rxn_params.alpha = 0.5;


[sol_total, params, rxn_params] = equation_solver(params, rxn_params);




function [value,isterminal,direction] = bounceEvents(t,y,yp,params)
%end simulation when certain mures is reached
%y(1)./(1-y(params.c_s_ind+1))-params_nondim.cap_max;
value = [y(end)-params.lower_volt;
    y(end)-params.upper_volt;
    sum(params.vfrac.*y(1:params.NM)) - params.cap_res_high;
    sum(params.vfrac.*y(1:params.NM)) - params.cap_res_low;];
% first one detects cutoff voltage, second one is the filling at the 
isterminal = [1; 1; 1; 1];   % Stop the integration
direction = [0; 0; 0; 0];   % Negative direction only
end



function centered = avg_2D(K, N, M)
%this function finds the centered values of a (c,T) matrix or a (T,c)
%matrix by reshaping, centering, and then re-outputing.
K_reshape = reshape(K, [N, M]);
%set side and boundary values to 0? depending on if it is c or T
centered = (K_reshape(1:end-1,:) + K_reshape(2:end,:))/2;
centered = centered(:);
end

function K_reshape = diff_2D(K, N, M, dx)
%this function finds the finite difference between values of a (c,T) matrix or a (T,c)
%matrix by reshaping, centering, and then re-outputing.
K_reshape = reshape(K, [N, M]);
%reshape matrix since ther are padded 0's along c or T
K_reshape = diff(K_reshape, 1, 1)/dx;
%K_reshape = diff(K_reshape, 1, 2)/dx;
K_reshape = K_reshape(:);
end


function y_reshape = pad_2D(y, N, M, method, value1, value2)
%pads function and reshapes so that calculations with Rxn and other terms
%are easier and value is the padding value 
y_reshape = reshape(y, N, M);
%pads vertical
if method == "R"
    %only pads c direction
    y_reshape = [zeros(N,1),y_reshape,zeros(N,1)];
elseif method == "R_left"
    %only pads c direction
    y_reshape = [zeros(N,1),y_reshape];
elseif method == "R_right"
    %only pads c direction
    y_reshape = [y_reshape,zeros(N,1)];
elseif method == "R_central"
    y_reshape = [value1,y_reshape,value2];
elseif method == "c_central"
    y_reshape = [value1.';y_reshape;value2.'];
elseif method == "c"
    y_reshape = [zeros(1,M);y_reshape;zeros(1,M)];
elseif method == "c_left"
    y_reshape = [zeros(1,M);y_reshape];
elseif method == "c_right"
    y_reshape = [y_reshape;zeros(1,M)];
end
y_reshape = y_reshape(:);
end


function gauss = gaussian(c_grid, R_grid, mu_c, sigma_c, mu_R, sigma_R)
%Returns points of a multivariate Gaussian in c_grid, T_grid space.
%c_grid, T_grid are the grid points of concentration and temperature, mu is
%the array of means of size 2, and sigma is the array of array of variances of size 2
%pos is an array constructed by packing the meshed arrays of variables
%x_1, x_2, x_3, ..., x_k into its _last_ dimension.
gauss = exp(-0.5*((c_grid-mu_c)/sigma_c).^2 -0.5*...
    ((R_grid-mu_R)/sigma_R).^2)/(sigma_c*sigma_R*2*pi);
end 



function out = mu_c(c, params)
%here we have the scaled chemical potential summed from the two different
%phases
%OCV from Colclasure 2020 NMC532
OCV = Colclasure_OCV(c);
%the OCV is in units of V, convert to kBT
out = -params.e/params.kBT*OCV + params.muR_ref;
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


function out = solve_mures(mures, IC, params, rxn_params, M)
%finds the inverse value that gives the current contraint we want, with
%mures = 0, which we autoset because of relaxation
%set T = 1
 %then we are solving for constant current
dfdt = dfdt_voltage(0, IC, mures, params, rxn_params);
out = M(1+params.NM+params.M,1:params.NM+params.M)*dfdt - params.control_value;
end


function out = logsig(x)
out = 1 ./ (1 + exp(-x));
end



function [k_array, rxn, eta] = R(c_grid, mures, k0, rxn_method, params)
%R is dc/dt, from the population balance equation. can be BV or linear
%set symmetry coefs
alpha = 0.5;
%get overpotential
muh = mu_c(c_grid, params).';
eta = muh-mures;
%for simplicity in calculation
c = c_grid.';
etaf = eta - log(c);
k_array = k0;
switch rxn_method
    case "linear"
        rxn = -k_array.*eta;
    case "BV"
        %using the transition state rxn model
        rxn = k_array.*(1-c).*(c./(1-c)).^alpha.*(exp(-alpha*eta)-...
            exp((1-alpha)*eta));
        if isa(eta,'casadi.MX') || isa(eta,'casadi.SX')
            eta = if_else(eta==0,1e-16,eta);
        else
            eta(eta==0) = 1e-16;
        end
        
        k_array = rxn./(-eta);
    case "Marcus"
        lmbda = params.lambda;
        %Marcus model
        rxn = k_array.*(1-c).*((exp(-power(lmbda+etaf,2)./(4*lmbda))) ...
          - (c.*exp(-power(lmbda-etaf,2)./(4*lmbda))));
        if isa(eta,'casadi.MX') || isa(eta,'casadi.SX')
            eta = if_else(eta==0,1e-16,eta);
        else
            eta(eta==0) = 1e-16;
        end
        k_array = rxn./(-eta);
    case "CIET"
        %temperoary
        i_red = helper_fun(-etaf, params.lambda);
        i_ox = helper_fun(etaf, params.lambda);
        rxn = k_array.*(1-c).*(i_red - c.*i_ox);
        if isa(eta,'casadi.MX') || isa(eta,'casadi.SX')
            eta = if_else(eta==0,1e-16,eta);
        else
            eta(eta==0) = 1e-16;
        end
        k_array = rxn./(-eta);
end
end



function rxn = R_new(c_grid, mures, v, k0, rxn_method, params)
%R is dc/dt, from the population balance equation. can be BV or linear
%set symmetry coefs
alpha = 0.5;
%get overpotential
muh = mu_c(c_grid./v, params).';
eta = muh-mures;
%for simplicity in calculation
c = c_grid.';
etaf = eta - log(c);
k_array = k0;
switch rxn_method
    case "linear"
        rxn = -k_array.*eta;
    case "BV"
        %using the transition state rxn model
        rxn = k_array.*(v.'-c).*(c./(v.'-c)).^alpha.*(exp(-alpha*eta)-...
            exp((1-alpha)*eta));
    case "Marcus"
        lmbda = params.lambda;
        %Marcus model
        rxn = k_array.*(v.'-c).*((exp(-power(lmbda+etaf,2)./(4*lmbda))) ...
          - (c.*exp(-power(lmbda-etaf,2)./(4*lmbda))));
    case "CIET"
        %temperoary
        i_red = helper_fun(-etaf, params.lambda);
        i_ox = helper_fun(etaf, params.lambda);
        rxn = k_array.*(v.'-c).*(i_red - c.*i_ox);
end
end





function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end


function out = dhelper_fundetaf(eta_f, lmbda)
out = (eta_f.*exp(-(lmbda - (eta_f.^2 + lmbda^(1/2) + 1).^(1/2)).^2./(4*lmbda)))...
    ./((exp(-eta_f) + 1).*(eta_f.^2 + lmbda^(1/2) + 1).^(1/2)) - ...
    (lmbda^(1/2)*pi^(1/2)*exp(-eta_f).*(erf((lmbda - ...
    (eta_f.^2 + lmbda^(1/2) + 1).^(1/2))./(2*lmbda^(1/2))) - 1))./(exp(-eta_f) + 1).^2;
end



function out = dideta(c_grid, mures, k0, rxn_method, params)
alpha = 0.5;
c = c_grid.';
muh = mu_c(c_grid, params).';
eta = muh-mures;
etaf = eta - log(c);
lmbda = params.lambda;
switch rxn_method
    case "BV"
        out = - k0.*(1-c).*(c./(1-c)).^alpha.*(-alpha*exp(-alpha*eta)-...
            (1-alpha)*exp((1-alpha)*eta));
    case "Marcus"
        %Marcus model
        out = k0.*((exp(-(etaf + lmbda).^2./(4.*lmbda)).*(2.*etaf + 2.*lmbda)) ...
            ./(4.*lmbda) - (c.*exp(-(etaf - lmbda).^2./(4*lmbda)).*(2.*etaf ...
            - 2.*lmbda))./(4.*lmbda)).*(c - 1);
    case "CIET"
        %temporary
        out = k0.*(1-c).*...
            (-dhelper_fundetaf(-etaf, lmbda) - c.*dhelper_fundetaf(etaf, lmbda));
end
end


function rxn = h(c_grid, mures, k0, rxn_method, params)
%R is dc/dt, from the population balance equation. can be BV or linear
%set symmetry coefs
alpha = 0.5;
%get overpotential
muh = mu_c(c_grid, params).';
eta = muh-mures;
%for simplicity in calculation
c = c_grid.';
etaf = eta - log(c);
k_array = k0;
switch rxn_method
    case "linear"
        rxn = -k_array.*eta;
    case "BV"
        %using the transition state rxn model
        rxn = k_array.*(exp(-alpha*eta)-...
            exp((1-alpha)*eta));
    case "Marcus"
        lmbda = params.lambda;
        %Marcus model
        rxn = k_array.*((exp(-power(lmbda+etaf,2)./(4*lmbda))) ...
          - (c.*exp(-power(lmbda-etaf,2)./(4*lmbda))));
    case "CIET"
        %temperoary
        i_red = helper_fun(-etaf, params.lambda);
        i_ox = helper_fun(etaf, params.lambda);
        rxn = k_array.*(i_red - c.*i_ox);
end
end


function out = dhdeta(c_grid, mures, k0, rxn_method, params)
alpha = 0.5;
c = c_grid.';
muh = mu_c(c_grid, params).';
eta = muh-mures;
etaf = eta - log(c);
lmbda = params.lambda;
switch rxn_method
    case "BV"
        out = - k0.*(-alpha*exp(-alpha*eta)-...
            (1-alpha)*exp((1-alpha)*eta));
    case "Marcus"
        %Marcus model
        out = - k0.*((exp(-(etaf + lmbda).^2./(4.*lmbda)).*(2.*etaf + 2.*lmbda)) ...
            ./(4.*lmbda) - (c.*exp(-(etaf - lmbda).^2./(4*lmbda)).*(2.*etaf ...
            - 2.*lmbda))./(4.*lmbda));
    case "CIET"
        %temporary
        out = k0.*...
            (-dhelper_fundetaf(-etaf, lmbda) - c.*dhelper_fundetaf(etaf, lmbda));
end
end




function out = Rinv(mures, c, control_value, k0, rxn_method, params)
[~, rxn, ~] = R(c,mures,k0,rxn_method,params);
out = abs(rxn-control_value);
end


function out = sum_2D(x_array, y_array, params)
out = zeros(params.M,1);
for i = 1:params.M
    %out(i) = trapz(x_array((i-1)*params.N+1:i*params.N), y_array((i-1)*params.N+1:i*params.N));
    out(i) = params.delta_c/params.N*(sum(y_array((i-1)*params.N+1:i*params.N))...
        -0.5*y_array((i-1)*params.N+1)-0.5*y_array(i*params.N));
end
end


function out = coth(x)
out = (exp(2*x)-1)./(exp(2*x)+1);
end


function out = Wfunc(params, Rxn_array, eta, omega, c_grid, mures, k0, rxn_method)
if params.order == 0
    out = 1./params.R_grid_1;
elseif params.order == 1
    out = (1+0.5*Rxn_array.'.*coth(0.5*eta.').*omega)./params.R_grid_1;
elseif params.order == 2
    out = (1+0.5*Rxn_array.'.*coth(0.5*eta.').*omega + ...
        (0.5*Rxn_array.').^2/2.*omega.^2)./params.R_grid_1;
elseif params.order == 3
    out = (1+0.5*Rxn_array.'.*coth(0.5*eta.').*omega + ...
        (0.5*Rxn_array.').^2/2.*omega.^2 + ...
        (0.5*Rxn_array.').^3/6.*coth(0.5*eta.').*omega.^3 ...
        )./params.R_grid_1; 
elseif params.order == 5
    out = (1+0.5*Rxn_array.'.*coth(0.5*eta.').*omega + ...
        (0.5*Rxn_array.').^2/2.*omega.^2 + ...
        (0.5*Rxn_array.').^3/6.*coth(0.5*eta.').*omega.^3 ...
        +(0.5*Rxn_array.').^4/24.*omega.^4 ... 
        +(0.5*Rxn_array.').^5/120.*coth(0.5*eta.').*omega.^5)./params.R_grid_1; 
elseif params.order == 10
    %this is MHC
    out = 1./params.R_grid_1.*((omega-c_grid)./(1-c_grid)).^0.5.*(1+...
        1./Rxn_array.'.*(c_grid.*(1-omega)./omega.^2).*dhdeta(c_grid, mures, k0, ...
        rxn_method, params).'.*dmudc(c_grid, params));
elseif params.order == -1
    %this is MHC
    out = 1./params.R_grid_1*1./(1-omega.*dideta(c_grid, mures, k0, ...
        rxn_method, params).');
elseif params.order == -2
    %this is MHC
    di = dideta(c_grid, mures, k0, ...
        rxn_method, params).';
    out = 1./params.R_grid_1*1./(1-omega.*di - omega.^2.*0.5.*(di.^2 + d2ideta2(c_grid, mures, k0, ...
        rxn_method, params).'));
elseif params.order == -10
    %this is MHC
    out = 1./params.R_grid_1.*((omega-c_grid)./(1-c_grid)).*(1+...
        1./Rxn_array.'.*(c_grid.*(1-omega)./omega.^2).*dhdeta(c_grid, mures, k0, ...
        rxn_method, params).'.*dmudc(c_grid, params));
end
end


function out = dmudc(c, params)
%here we have the scaled chemical potential summed from the two different
%phases
%OCV from Colclasure 2020 NMC532
OCV = dColclasure_OCVdx(c);
%the OCV is in units of V, convert to kBT
out = -params.e/params.kBT*OCV + params.muR_ref;
end



function OCV = dColclasure_OCVdx(x)
%colclasure OCV function
OCV = -3.640117692001490E+03*14*x.^13.0 + 1.317657544484270E+04*13*x.^12.0...
- 1.455742062291360E+04*12*x.^11.0 - 1.571094264365090E+03*11*x.^10.0...
+ 1.265630978512400E+04*10*x.^9.0 - 2.057808873526350E+03*9*x.^8.0...
- 1.074374333186190E+04*8*x.^7.0 + 8.698112755348720E+03*7*x.^6.0...
- 8.297904604107030E+02*6*x.^5.0 - 2.073765547574810E+03*5*x.^4.0...
+ 1.190223421193310E+03*4*x.^3.0 - 2.724851668445780E+02*3*x.^2.0...
+ 2.723409218042130E+01*2*x - 4.158276603609060E+00 +...
-5.573191762723310E-04*exp(6.560240842659690E+00*x.^4.148209275061330E+01)*6.560240842659690E+00*4.148209275061330E+01.*x.^3.148209275061330E+01;
end


function out = d2helper_fundetaf2(eta_f, lmbda)
out = exp(-(lmbda - (eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)).^2./(4*lmbda))./...
    ((exp(-eta_f) + 1).*(eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)) + ...
    (3991211251234741*lmbda.^(1/2).*exp(-eta_f).*(erf((lmbda - (eta_f.^2 + ...
    lmbda.^(1/2) + 1).^(1/2))./(2.*lmbda^(1/2))) - 1))./(2251799813685248.*(exp(-eta_f) + 1).^2)...
    - (3991211251234741*lmbda^(1/2).*exp(-2*eta_f).*(erf((lmbda - (eta_f.^2 ...
    + lmbda^(1/2) + 1).^(1/2))./(2*lmbda^(1/2))) - 1))./(1125899906842624.*...
    (exp(-eta_f) + 1).^3) - (eta_f.^2.*exp(-(lmbda - (eta_f.^2 + lmbda^(1/2)...
    + 1).^(1/2)).^2./(4*lmbda)))./((exp(-eta_f) + 1).*(eta_f.^2 + lmbda^(1/2)...
    + 1).^(3/2)) + (eta_f.*exp(-eta_f).*exp(-(lmbda - (eta_f.^2 + lmbda^(1/2)...
    + 1).^(1/2)).^2./(4*lmbda)))./((exp(-eta_f) + 1).^2.*(eta_f.^2 + lmbda^(1/2)...
    + 1).^(1/2)) + (eta_f.^2.*exp(-(lmbda - (eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)).^2. ...
    ./(4*lmbda)).*(lmbda - (eta_f.^2 + lmbda^(1/2) + 1).^(1/2)))./(2*lmbda.*(exp(-eta_f)...
    + 1).*(eta_f.^2 + lmbda^(1/2) + 1)) + (3991211251234741.*eta_f.*exp(-eta_f) ...
    .*exp(-(lmbda - (eta_f.^2 + lmbda^(1/2) + 1).^(1/2)).^2./(4*lmbda)))./ ...
    (2251799813685248*pi^(1/2).*(exp(-eta_f) + 1).^2.*(eta_f.^2 + lmbda^(1/2) + 1).^(1/2));
end



function out = d2ideta2(c_grid, mures, k0, rxn_method, params)
alpha = 0.5;
c = c_grid.';
muh = mu_c(c_grid, params).';
eta = muh-mures;
etaf = eta - log(c);
lmbda = params.lambda;
switch rxn_method
    case "BV"
        out = -k0.*(exp(-eta.*(alpha - 1)).*(alpha - 1).^2 - alpha.^2 ...
        .*exp(-alpha.*eta)).*(c - 1).*(-c./(c - 1)).^alpha;
    case "Marcus"
        %Marcus model
        out = k0.*(c - 1).*(exp(-(etaf + lmbda).^2./(4*lmbda))./(2.*lmbda) ...
            - (c.*exp(-(etaf - lmbda).^2./(4*lmbda)))./(2*lmbda) - ...
            (exp(-(etaf + lmbda).^2./(4.*lmbda)).*(etaf + lmbda).^2)./(4*lmbda.^2) ...
            + (c.*exp(-(etaf - lmbda).^2./(4*lmbda)).*(etaf - lmbda).^2)./(4*lmbda^2));
    case "CIET"
        %temporary
        out = k0.*(1-c).*...
            (d2helper_fundetaf2(-etaf, lmbda) - c.*d2helper_fundetaf2(etaf, lmbda));
end
end






function [df_dt, y_array] = dfdt_voltage(t, y, mures, params, rxn_params)
%Generate dfdt for each steps, using R, dTdt with 0 padding, the c and T discretizations
%and the diffusive terms for concentration and temperature, D at CENTERED values, and kappa
%returns df_dt
%driving force of concentration at the center of each point
%f*R term, padded with zeros since taken only at center value

%to find the midpoints between the values, we average values
f = y(1:params.NM);
v = y(params.NM+1:params.NM+params.M);
%reshape and stack into NM from
v = reshape(repmat(v,1,params.N)',params.NM,1);

% first, test for the terms that are blowing up
[k_array, Rxn_array, eta] = R(params.c_grid_1, mures, rxn_params.k0, ...
    params.rxn_method, params);
Rxn_total = R_new(params.c_grid_1, mures, v, rxn_params.k0, ...
    params.rxn_method, params);
% if isa(Rxn_total,'casadi.MX') || isa(Rxn_total,'casadi.SX')
%     Rxn_total = if_else(abs(Rxn_total) < 1e-3<0,sign(Rxn_total),Rxn_total);
% else
%     Rxn_total(abs(Rxn_total) < 1e-3) = sign(Rxn_total(abs(Rxn_total) < 1e-3));
% end
% if isa(Rxn_array,'casadi.MX') || isa(Rxn_array,'casadi.SX')
%     Rxn_array = if_else(abs(Rxn_array) < 1e-3<0,sign(Rxn_array),Rxn_array);
% else
%     Rxn_array(abs(Rxn_array) < 1e-3) = sign(Rxn_array(abs(Rxn_array) < 1e-3));
% end

% y_array = R_resistance(mures, Rxn_array.', rxn_params.k0_res, ...
%     params).*f.*params.c_grid_1;
y_array = R_resistance(mures, Rxn_total.', rxn_params.k0_res, ...
    params).*f.*params.c_grid_1./params.R_grid_1;
domegadt = [];
for i = 1:params.M
    domegadt = [
        domegadt
        -(params.delta_c*y_array((i-1)*params.N+1)/2 + params.delta_c*...
        sum(y_array((i-1)*params.N+2:(i-1)*params.N+params.N-1)) + ...
        params.delta_c*y_array(params.N*i)/2);
        ];
end

    
domegadt = domegadt./params.IC1DR.';

%remove the indices where the system blows up

%next, get the diffusion equations
% W = Wfunc(params, Rxn_array, eta, v, params.c_grid_1, mures, rxn_params.k0, params.rxn_method);
% %remove the indexes that are insanely large
% if isa(W,'casadi.MX') || isa(W,'casadi.SX')
%     W = if_else(W<0,0,W);
% else
%     W(W<0) = 0;
% end
%for numerical efficiency, where f is very small, we can drop the W terms
F_c = avg_2D(f.*Rxn_total.'./params.R_grid_1,params.N,params.M);
% F_c = avg_2D(f.*Rxn_array.'.*W,params.N,params.M);
%pad with zeros and add diffusive term
F_c = F_c - rxn_params.D0*avg_2D(k_array,params.N,params.M).* ...
    diff_2D(f./params.R_grid_1,params.N,params.M,params.delta_c);
% F_c = F_c - rxn_params.D0*avg_2D(k_array,params.N,params.M).* ...
%     diff_2D(f.*W,params.N,params.M,params.delta_c);
F_c = pad_2D(F_c,params.N-1,params.M,"c");

%set final equation
% df_dt(1:params.NM,1) = -diff_2D(F_c,params.N+1,params.M,params.delta_c);
df_dt = -diff_2D(F_c,params.N+1,params.M,params.delta_c);

df_dt = [df_dt; domegadt];


%         (sum([params.delta_c/2; ones(params.N-2,1); params.delta_c/2].*y_array((i-1)*params.N+1:params.N*i)));
   
%         (sum([params.delta_c/2; ones(params.N-2,1); params.delta_c/2].*f((i-1)*params.N+1:params.N*i).*y_array((i-1)*params.N+1:params.N*i)))...
%         ./(sum([params.delta_c/2; ones(params.N-2,1); params.delta_c/2].*f((i-1)*params.N+1:params.N*i)));
%     

%         (sum(y_array((i-1)*params.N+1:i*params.N))-0.5*y_array((i-1)*params.N+1)-0.5*y_array(i*params.N))...
%         /(sum(f((i-1)*params.N+1:i*params.N))-0.5*f((i-1)*params.N+1)-0.5*f(i*params.N));
% %since this growth is per area.
end


function df_dt = dfdt_current(t, y, params, rxn_params)
%Generate dfdt for each steps, using R, dTdt with 0 padding, the c and T discretizations
%and the diffusive terms for concentration and temperature, D at CENTERED values, and kappa
%returns df_dt
%driving force of concentration at the center of each point
%f*R term, padded with zeros since taken only at center value
%split y values
mures = y(end);
y = y(1:end-1);

%use constant potential wrapper
[df_dt, y_array] = dfdt_voltage(t, y, mures, params, rxn_params);

% add current constraint
%df_dt(end+1,1) = params.control_value - sum(params.vfrac.*y_array./params.R_grid_1);
df_dt = [df_dt; params.control_value];
end


function rxn = R_resistance(mures, R_intercalation, k, params)
%R is dc/dt, from the population balance equation
%eta = params.V_res - mures - R_intercalation.*Omega;
eta = params.V_res - mures;
%BV
k_array = ones(size(eta))*k;
%using the transition state rxn model
rxn = k_array.*(exp(eta));
end


function gauss = gaussian1D(c_grid, mu, sigma)
%Returns points of a multivariate Gaussian in c_grid, T_grid space.
%c_grid, T_grid are the grid points of concentration and temperature, mu is
%the array of means of size 2, and sigma is the array of array of variances of size 2
%pos is an array constructed by packing the meshed arrays of variables
%x_1, x_2, x_3, ..., x_k into its _last_ dimension.
gauss = exp(-0.5*((c_grid-mu)/sigma).^2)/(sigma*sqrt(2*pi));
end 


function [sol_total, params, rxn_params] = equation_solver(params, rxn_params)
%runs ODE: with inputs final timestep, c and T discretization, temperature
%range, and mures we are applying
%set parameters

%takes in the final timestep, the c discretization, and the set chemical
%potential
%returns the discretizations of c an time, as well as teh solution in (c,t)

%set initial conditions
%get initial_condition with mu, sigma
IC1Dc = gaussian1D(params.c_array, params.c0, params.sigma_c);
sumc = trapz(params.delta_c, IC1Dc);
IC1Dc = IC1Dc/sumc;
%renormalize this
IC1DR = gaussian1D(params.R_array, params.avg_R, params.sigma_R); 
%for numerical efficiency, we don't renormalize IC1DR until later.
% sumR = trapz(params.delta_R, IC1DR);
% IC1DR = IC1DR/sumR;
sumR = 1e6;
%sumR = 1e9;
IC1DR = IC1DR/sumR;
IC = gaussian(params.c_grid, params.R_grid, params.c0, params.sigma_c, params.avg_R, params.sigma_R);
%IC = IC/(sumR*sumc);
%IC = IC/sumc;
IC = IC/sumR;

params.IC1DR = IC1DR;

%first stack in c arrays. then stack in R arrays.

NM = params.N*params.M;
N = params.N;
M = params.M;
M_mat = spdiags(ones(NM+M,1),[0],NM+M,NM+M);
% av_vols = sum(params.R_array.^3)./params.M;
% %padding for the edges being 1/2 the value of the central values
% params.vfrac = pad_2D(ones(N-2,M)*params.delta_c,N-2,M,'c_central',...
%     params.delta_c/2*ones(M,1),params.delta_c/2*ones(M,1)).*...
%     params.R_grid_1.^3./av_vols;
total_vol = trapz(params.delta_R, params.R_array.^3.*IC1DR);
%padding for the edges being 1/2 the value of the central values
params.vfrac = pad_2D(ones(N-2,M)*params.delta_c,N-2,M,'c_central',...
    params.delta_c/2*ones(M,1),params.delta_c/2*ones(M,1)).*...
    pad_2D(ones(N,M-2)*params.delta_R,N,M-2,'R_central',...
    params.delta_R/2*ones(N,1),params.delta_R/2*ones(N,1)).*...
    params.R_grid_1.^3./total_vol;

% %test
% IC = zeros(size(IC));
% IC(1:params.N) = IC1Dc;

%then turn into array
IC = IC(:);
%here, append the zero values for the initial values of omega
IC = [IC; ones(params.M,1)];

switch params.protocol
    case "voltage"
        %jacobian is constant

        %gets jacobian
%         disp("calculate Jacobian")
%         tic()
%         jac = jacobian_voltage(R_array,dTdt_array,k_array,T_grid_1,delta_c,...
%             delta_T,N,M,D0);
%         toc()
        
        %solve ODE keeping positive
        %opts = odeset('Jacobian',jac);
        opts = odeset('NonNegative',1:length(c_array)*length(T_array));
        disp("solve ODE")
        tic()
        sol = ode113(@(t,f) dfdt_voltage(t,f,params.control_value,params,rxn_params));
        %[t_array, f_solution] = ode15s(@(t,f) dfdt_voltage(t,f,R_array,...
        %    dTdt_array,k_array,T_grid_1,delta_c,delta_T,N,M,D0),...
        %    [0,final_t], IC, opts);
        toc()

        %get solution only at time points that we want
        f_solution = reshape(deval(sol,t_array),[num_times,params.N,params.M]);
        
        
    case "current"
        %here, y is longer by one value because last value is mures
        %append value to IC
        
        %get solution time points
        num_times = 50;

        %mass matrix for trapezoidal integration
        M_mat(NM+M+1,1:NM) = params.c_grid_1.*params.vfrac; % for the intercalation current
        %now we add in the degradation current
        M_mat(NM+M+1,NM+1:NM+M) = [params.delta_R/2, params.delta_R*ones(1,M-2), params.delta_R/2].*...
            params.R_array.^2.*IC1DR./params.beta./total_vol;
        
        %padding for the edges being 1/2 the value of the central values
%        M_mat(NM+M+1,NM+1:NM+M) = 1/(params.beta*params.R_array); % for the intercalation current
        M_mat(NM+M+1,NM+M+1) = 0;

        %find initial guess of mures
        options = optimset('Display','on','TolX',1e-10,'TolFun',1e-10);

        %we only charge up to 1 for constant current
        final_t = params.final_t;
        t_array = linspace(0,final_t,num_times);
             
        t_end = 0;
        sol_total = {};
        sol_total.x = [];
        sol_total.y = [];
        times_charge = [];
        times_discharge = [];
        
        muresguess1 = 0;
        muresguess2 = 20;
        
        if params.voltage_cutoffs
            opts = odeset('Mass', M_mat, 'Jacobian', @(t,u) ...
                AutoDiffJacobianAutoDiff(@(x) dfdt_current(t,x,params,rxn_params), u), ...
                'RelTol', 1e-9, 'AbsTol', 1e-9, 'MStateDependence', 'none');
        else
            opts = odeset('Mass', M_mat, 'Jacobian', @(t,u) ...
                AutoDiffJacobianAutoDiff(@(x) dfdt_current(t,x,params,rxn_params), u), ...
                'RelTol', 1e-9, 'AbsTol', 1e-9, 'Events', @(t,y) bounceEvents(t, y, params), ...
                'MStateDependence', 'none');
        end
        
        %finish setting mass matrix
        for i = 1:params.num_cycles
            i
% lets try to cycle 
            params.control_value = abs(params.control_value);
            mures0 = fmincon(@(mures) Rinv(mures, params.c0, params.control_value,...
                rxn_params.k0, params.rxn_method, params),muresguess1);
            %find mures that gives consistent IC
            %interpolate to find the value that is 0 or closest to 0
            %[mures,~,exitflag,~] = fsolve(@(mures) solve_mures(mures, IC, params, rxn_params, M_mat), muresguess1, options);
            [mures,~,exitflag,~] = fzero(@(mures) solve_mures(mures, IC, params, rxn_params, M_mat), mures0, options);
            % dfdt = dfdt_voltage(0, IC, mures, params, rxn_params);
            muresguess1 = mures;
            IC = [IC; mures];

            % create a fast version of the dfdt function and Jacobian in casadi
            addpath('~/codes/casadi')
            import casadi.*
            xsym    = SX.sym('x',[length(IC),1]);
            tsym    = SX.sym('t',1);
            % Get the symbolic model equations
            df_tot = dfdt_current(0,xsym,params,rxn_params);
        
            dfdt_current_tmp = Function('f',{tsym,xsym}, {df_tot});
            dfdt_current_fast = @(t,x) full(dfdt_current_tmp(t,x));
            % Evaluate the Jacobian matrix
            J = jacobian(df_tot, xsym);
            
            jacobian_tmp = Function('fJ',{tsym, xsym}, {J});
            jacobian_fast = @(t,x) sparse(jacobian_tmp(t,x));

            if ~params.voltage_cutoffs
                opts = odeset('Mass', M_mat, 'Jacobian', jacobian_fast, ...
                    'RelTol', 1e-8, 'AbsTol', 1e-8);
            else
                opts = odeset('Mass', M_mat, 'Jacobian', jacobian_fast, ...
                    'RelTol', 1e-7, 'AbsTol', 1e-8, 'Events', @(t,y) bounceEvents(t, y, params));
            end   
           
           % tic()
           % sol = ode15s(@(t,f) dfdt_current_fast(t,f), [0,final_t], IC, opts);
           % toc()

            % implicit
             if ~params.voltage_cutoffs
                 opts = odeset('Mass', M_mat, 'MassSingular', 'yes', 'Jacobian', @(t,u,du) jac_implicit(t,u,du,jacobian_fast,M_mat), ...
                     'RelTol', 1e-8, 'AbsTol', 1e-8);
             else
                 opts = odeset('Mass', M_mat, 'MassSingular', 'yes', 'Jacobian', @(t,u,du) jac_implicit(t,u,du,jacobian_fast,M_mat), ...
                     'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y,yp) bounceEvents(t, y, yp, params));
             end
 
             df=dfdt_current_fast(0,IC);
             tic()
            % final_t = 5e-3;
             sol = ode15i(@(t,f,df) dfdt_current_fast(t,f) - M_mat*df, [0,final_t], IC, df, opts);
             toc()
             IC = [];
%             tic()
%             sol = ode15s(@(t,f) dfdt_current(t,f,params,rxn_params), [0,final_t], IC, opts);
%             toc()

            times_charge = [times_charge, sol.x(end)];
            t_array = linspace(0,sol.x(end),num_times);
            y = deval(sol,t_array);
            sol_total.x = [sol_total.x, t_array+t_end];
            sol_total.y = [sol_total.y, y];
            IC = sol_total.y(1:end-1,end);
            t_end = sol_total.x(end);
            
            params.control_value = -abs(params.control_value);
            % need a new function since the parameters changed
            % Get the symbolic model equations
            if ~params.voltage_cutoffs
                opts = odeset('Mass', M_mat, 'MassSingular', 'yes', 'Jacobian', @(t,u,du) jac_implicit(t,u,du,jacobian_fast,M_mat), ...
                    'RelTol', 1e-8, 'AbsTol', 1e-8);
            else
                opts = odeset('Mass', M_mat, 'MassSingular', 'yes', 'Jacobian', @(t,u,du) jac_implicit(t,u,du,jacobian_fast,M_mat), ...
                    'RelTol', 1e-8, 'AbsTol', 1e-8, 'Events', @(t,y,yp) bounceEvents(t, y, yp, params));
            end

            mures0 = fmincon(@(mures) Rinv(mures, 0.9, params.control_value,...
               rxn_params.k0, params.rxn_method, params), muresguess2);
            %find mures that gives consistent IC
            %interpolate to find the value that is 0 or closest to 0
            
            %[mures,~,exitflag,~] = fsolve(@(mures) solve_mures(mures, IC, params, rxn_params, M_mat), muresguess2, options);
            [mures,~,exitflag,~] = fzero(@(mures) solve_mures(mures, IC, params, rxn_params, M_mat), mures0, options);
            muresguess2 = mures;
            IC = [IC; mures];

            xsym    = SX.sym('x',[length(IC),1]);
            tsym    = SX.sym('t',1);
            % Get the symbolic model equations
            df_tot = dfdt_current(0,xsym,params,rxn_params);
        
            dfdt_current_tmp = Function('f',{tsym,xsym}, {df_tot});
            dfdt_current_fast = @(t,x) full(dfdt_current_tmp(t,x));
            % Evaluate the Jacobian matrix
            J = jacobian(df_tot, xsym);
            
            jacobian_tmp = Function('fJ',{tsym, xsym}, {J});
            jacobian_fast = @(t,x) sparse(jacobian_tmp(t,x));

            df=dfdt_current_fast(0,IC);
            tic()
            sol = ode15i(@(t,f,df) dfdt_current_fast(t,f) - M_mat*df, [0,final_t], IC, df, opts);
            toc()
            IC = [];
%             tic()
%             sol = ode15s(@(t,f) dfdt_current_fast(t,f), [0,final_t], IC, opts);
%             toc()
%             tic()
%             sol = ode15s(@(t,f) dfdt_current(t,f,params,rxn_params), [0,final_t], IC, opts);
%             toc()
            times_discharge = [times_discharge, sol.x(end)];
            t_array = linspace(0,sol.x(end),num_times);
            y = deval(sol,t_array);
            sol_total.x = [sol_total.x, t_array+t_end];
            sol_total.y = [sol_total.y, y];
            IC = sol_total.y(1:end-1,end);
            t_end = sol_total.x(end);
            
            save('fitness_distribution.mat', 'sol_total', 'params', 'rxn_params', 'times_discharge', 'times_charge', '-v7.3')
        end
       
end
% 
%  figures
% 
%figure 1: charge/discharge effect on radial distribution effect
%plot solution
% 
% 
%s

%save('fitness_distribution.mat', 'sol_total', 'params', 'rxn_params', 'times_discharge', 'times_charge')


subplot(2,1,1)

title('Concentration')
%c = get_total_concentration(sol_total.x, sol_total.y(1:params.NM,:), params);
 [m, v] = get_mean_variance(sol_total.y(1:params.NM,:), params);

image([0,sol_total.x(end)], [min(params.R_array), max(params.R_array)]*100, m, 'CDataMapping', 'scaled')
colorbar
%here are the labels we want for average concentration
t_ticks = [0, final_t/2, final_t, 3*final_t/2, 2*final_t];
c_ticks = [0.4, 0.65, 0.9, 0.65, 0.4];
set(gca,'xtick',t_ticks,'xticklabel',c_ticks)
xlabel('State of Charge')
ylabel('Particle Radius (nm)')
subplot(2,1,2)
title('Standard Deviation of Concentration')
image([0,sol_total.x(end)], [min(params.R_array), max(params.R_array)]*100, v, 'CDataMapping', 'scaled')
colorbar
set(gca,'xtick',t_ticks,'xticklabel',c_ticks)
xlabel('State of Charge')
ylabel('Particle Radius (nm)')

savefig('fig1.fig')

figure(2)

% smooth data
newpoints = 100;
[xq,yq] = meshgrid(...
            linspace(0, params.num_cycles, newpoints ),...
            linspace(min(params.R_array)*100, max(params.R_array)*100, newpoints)...
          );
%get the unique values
[xval, ia, ic] = unique(sol_total.x/max(sol_total.x) * (params.num_cycles));
yval = sol_total.y(params.NM+1:params.NM+params.M,:)./min(sol_total.y(params.NM+1:params.NM+params.M,:));
yval = yval(:,ia);

BDmatrixq = interp2(xval.', params.R_array*100, yval ,xq,yq,'cubic');
[c,h]=contourf(xq,yq,BDmatrixq);
set(h,'LineColor','none');

% image([0, params.num_cycles], [min(params.R_array), max(params.R_array)]*100, ...
%     sol_total.y(params.NM+1:params.NM+params.M,:)./min(sol_total.y(params.NM+1:params.NM+params.M,:)), 'CDataMapping', 'scaled')
xlabel('Cycle Number')
ylabel('Particle Radius (nm)')
colorbar
savefig('fig2.fig')


figure(3)
%average amount of time spent at "low potential regions"
[m, v] = get_mean_variance(sol_total.y(1:params.NM,:), params);
hits = get_low_hits(sol_total.y(1:params.NM,:), m, params);
image([0,params.num_cycles], [min(params.R_array), max(params.R_array)]*100, cumsum(hits,2), 'CDataMapping', 'scaled');
xlabel('Cycle Number')
ylabel('Particle Radius (nm)')
colorbar
savefig('fig3.fig')

figure(4)
subplot(2,1,1)
newpoints = 500;
%PLOT INITIal and final fitness landscapes
[k_array, Rxn_array, eta] = R(params.c_grid_1, sol_total.y(end,1), rxn_params.k0, params.rxn_method, params);
omega = sol_total.y(params.NM+1:params.NM+params.M,1);
W_init = reshape(Wfunc(params, Rxn_array, eta, reshape(repmat(omega,1,params.N)',params.NM,1), params.c_grid_1, sol_total.y(end,1), rxn_params.k0, params.rxn_method), params.N, params.M);
%W_init = reshape(Wfunc(params, Rxn_array, eta, reshape(repmat(omega,1,params.N)',params.NM,1)), params.N, params.M);
[k_array, Rxn_array, eta] = R(params.c_grid_1, sol_total.y(end,end), rxn_params.k0, params.rxn_method, params);
omega = sol_total.y(params.NM+1:params.NM+params.M,end);
W_final = reshape(Wfunc(params, Rxn_array, eta, reshape(repmat(omega,1,params.N)',params.NM,1), params.c_grid_1, sol_total.y(end,end), rxn_params.k0, params.rxn_method), params.N, params.M);
%W_final = reshape(Wfunc(params, Rxn_array, eta, reshape(repmat(omega,1,params.N)',params.NM,1)), params.N, params.M);
% smooth data
%image([min(params.R_array), max(params.R_array)]*100, [min(params.c_array), max(params.c_array)], W_init, 'CDataMapping', 'scaled')
[xq,yq] = meshgrid(...
            linspace(min(params.R_array)*100, max(params.R_array)*100, newpoints),...
            linspace(min(params.c_array), max(params.c_array), newpoints),...
            newpoints);
%get the unique values
zq = interp2(params.R_array*100, params.c_array, W_init./W_init(:,end), xq, yq, 'cubic');
[c,h]=contourf(xq,yq,zq);
set(h,'LineColor','none');
ylabel('Concentration')
xlabel('Particle Radius (nm)')
colorbar
set(gca, 'clim', [0, 30]);

% image([0, params.num_cycles], [min(params.R_array), max(params.R_array)]*100, ...
%     sol_total.y(params.NM+1:params.NM+params.M,:)./min(sol_total.y(params.NM+1:params.NM+params.M,:)), 'CDataMapping', 'scaled')
subplot(2,1,2)
%image([min(params.R_array), max(params.R_array)]*100, [min(params.c_array), max(params.c_array)], W_final, 'CDataMapping', 'scaled')
zq = interp2(params.R_array*100, params.c_array, W_final./W_final(:,end), xq, yq, 'cubic');
[c,h]=contourf(xq,yq,zq);
% image([0, params.num_cycles], [min(params.R_array), max(params.R_array)]*100, ...
%     sol_total.y(params.NM+1:params.NM+params.M,:)./min(sol_total.y(params.NM+1:params.NM+params.M,:)), 'CDataMapping', 'scaled')
ylabel('Concentration')
xlabel('Particle Radius (nm)')
colorbar
set(gca, 'clim', [0, 30]);
set(h,'LineColor','none');
savefig('fig4.fig')


figure(5)
%time spent in each cycle
plot(1:1:params.num_cycles, times_charge/times_charge(1))
hold on
plot(1:1:params.num_cycles, times_discharge/times_discharge(1))
legend('Charge', 'Discharge')
xlabel('Cycle Number')
ylabel('Nondimensionalized Time')
savefig('fig5.fig')


end

function [du,dup] = jac_implicit(t,u,du,jacobian_fast,M_mat)
    du = jacobian_fast(t,u);
    dup = -M_mat;
end


function av = get_total_concentration(t, f, params)
%returns total average concentration in the form of a polynomial
%(concentration=x, time=y)
av = (params.vfrac.*params.c_grid_1).'*f;
end


function [m, v] = get_mean_variance(f, params)
%return the mean, variance for each R value
f_prime = reshape(f, [params.N, params.M, size(f,2)]);
m = sum(f_prime.*[params.delta_c/2; ones(params.N-2,1); ...
    params.delta_c/2].*params.c_grid)./sum(f_prime.*[params.delta_c/2; ...
    ones(params.N-2,1); params.delta_c/2]);
v = squeeze(sum(f_prime.*[params.delta_c/2; ones(params.N-2,1); ...
    params.delta_c/2].*(params.c_grid-m).^2)./sum(f_prime.*[params.delta_c/2; ...
    ones(params.N-2,1); params.delta_c/2]));
m = squeeze(m);
end

function hits = get_low_hits(f, m, params)
%get the number of hits we spend at low time
f_prime = reshape(f, [params.N, params.M, size(f,2)]);
save_ind = repmat((params.c_grid < 0.38), 1, 1, size(f,3));
f_prime = f_prime.*save_ind;
hits = squeeze(sum(f_prime.*[params.delta_c/2; ones(params.N-2,1); ...
    params.delta_c/2].*params.c_grid)./sum(f_prime.*[params.delta_c/2; ...
    ones(params.N-2,1); params.delta_c/2]))./m;
end


