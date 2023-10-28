
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




figure(6) %iredoi fraction and 
%copy paste reaction rate and Wfunc. taken at 5, 10, 100 mV
% which is 0.2 kBT, 0.4 kBT, 4 kBT
eta = linspace(-8,8);
%taken at kBT = 5, 10, 100
plot(eta, dideta(0.3, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end,:))
hold on
plot(eta, dideta(0.5, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end-2,:))
plot(eta, dideta(0.7, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end-4,:))
xline(0.2, '--')
xline(-0.2, '--')
xline(0.8, '--')
xline(-0.8, '--')
xline(4, '--')
xline(-4, '--')

[x_red, y_red] = get_red_points([0.2, 0.8, 4], [0.3, 0.5, 0.7], 10, 3.78);
scatter(x_red, y_red, 36, 'red', 'filled')
% semilogy(eta, dideta(0.3, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end,:))
% hold on
% semilogy(eta, dideta(0.5, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end-2,:))
% semilogy(eta, dideta(0.7, eta, 10, 3.78), 'LineWidth',2, 'Color', map_blue_less(end-4,:))
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('c = 0.3', 'c = 0.5', 'c = 0.7', 'Fontsize', 20)
% legend boxoff
export_fig fig26.png -transparent -m3 -r300



figure(7) %dideta_values
%copy paste reaction rate and Wfunc
eta = linspace(-8,8);
plot(eta, iredoi(0.3, eta, 3.78),'LineWidth',2, 'Color', map_blue_less(end,:))
hold on
plot(eta, iredoi(0.5, eta, 3.78),'LineWidth',2, 'Color', map_blue_less(end-2,:))
plot(eta, iredoi(0.7, eta, 3.78),'LineWidth',2, 'Color', map_blue_less(end-4,:))
xline(0.2, '--')
xline(-0.2, '--')
xline(0.8, '--')
xline(-0.8, '--')
xline(4, '--')
xline(-4, '--')
[x_red, y_red] = get_red_points2([0.2, 0.8, 4], [0.3, 0.5, 0.7], 3.78);
scatter(x_red, y_red, 36, 'red', 'filled')
ylim([-2, 2])
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
%HPPC EXAMPLE
export_fig fig27.png -transparent -m3 -r300


figure(8) %i current
%copy paste reaction rate and Wfunc
semilogy(eta, abs(i(0.3, eta, 10, 3.78)), 'LineWidth',2, 'Color', map_blue_less(end,:))
hold on
semilogy(eta, abs(i(0.5, eta, 10, 3.78)), 'LineWidth',2, 'Color', map_blue_less(end-2,:))
semilogy(eta, abs(i(0.7, eta, 10, 3.78)), 'LineWidth',2, 'Color', map_blue_less(end-4,:))
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
legend('c = 0.3', 'c = 0.5', 'c = 0.7', 'Fontsize', 20, 'Location', 'southeast')
legend boxoff
% legend('c = 0.3', 'c = 0.5', 'c = 0.7', 'Fontsize', 20)
% legend boxoff
export_fig fig28.png -transparent -m3 -r300


figure(9) %i current
%copy paste reaction rate and Wfunc
ctilde = linspace(0.8,1);
plot(ctilde, (ctilde-0.3)/(1-0.3), 'LineWidth',2, 'Color', map_blue_less(end,:))
hold on
plot(ctilde, (ctilde-0.5)/(1-0.5), 'LineWidth',2, 'Color', map_blue_less(end-2,:))
plot(ctilde, (ctilde-0.7)/(1-0.7), 'LineWidth',2, 'Color', map_blue_less(end-4,:))
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('c = 0.3', 'c = 0.5', 'c = 0.7', 'Fontsize', 20, 'Location', 'southeast')
% legend boxoff
% legend('c = 0.3', 'c = 0.5', 'c = 0.7', 'Fontsize', 20)
% legend boxoff
export_fig fig29.png -transparent -m3 -r300





figure(1)
%copy paste reaction rate and Wfunc
semilogy(linspace(0,0.25,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_R_f_c_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(:,1),'LineWidth',2, 'Color', map_purple_less(end,:))
hold on
load('../full_cell_pulse_fit_R_f_c_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(:,1),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_R_f_c_only_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(:,1),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
xlim([1,10])
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
export_fig fig21.png -transparent -m3 -r300


figure(2)
%copy paste reaction rate and Wfunc
plot(linspace(1,0.9,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_c_tilde_c_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(:,2),'LineWidth',2, 'Color', map_purple_less(end-4,:))
hold on
load('../full_cell_pulse_fit_c_tilde_c_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(:,2),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_c_tilde_c_only_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(:,2),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
xlim([1,10])
export_fig fig22.png -transparent -m3 -r300


figure(3)
%copy paste reaction rate and Wfunc
plot(linspace(1000,900,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--', 'HandleVisibility','off')
hold on
load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(:,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'HandleVisibility','off')
hold on
load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(:,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'LineStyle', ':', 'HandleVisibility','off')
load('../current_pulse_results/full_cell_pulse_fit_c_lyte_only_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(:,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'LineStyle', '-.', 'HandleVisibility','off')
% load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_charge_linear.mat')
% plot(zPulse(:,3)*1000,'o','LineWidth',2, 'Color', map_purple_less(end-4,:))
% hold on
% load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_discharge_linear.mat')
% plot(zPulse(:,3)*1000,'x','LineWidth',2, 'Color', map_purple_less(end-4,:))
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
xlim([1,10])
% legend('Discharge Fitted Solution' , 'Charge Fitted Solution', 'Fontsize', 20, 'Location', 'southwest')
% legend boxoff
%HPPC EXAMPLE
export_fig fig23.png -transparent -m3 -r300


figure(4)
%copy paste reaction rate and Wfunc
semilogy(linspace(0,0.25,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_R_f_a_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(:,4),'LineWidth',2, 'Color', map_purple_less(end,:))
load('../full_cell_pulse_fit_R_f_a_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(:,4),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_R_f_a_only_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(:,4),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', '-.')

set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
xlim([1,10])
export_fig fig24.png -transparent -m3 -r300

f = figure(5)
%copy paste reaction rate and Wfunc
plot(linspace(1,0.9,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_c_tilde_a_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:))
load('../full_cell_pulse_fit_c_tilde_a_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_c_tilde_a_only_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
legend('Exact Solution', 'Voltage Pulse Linear Fit', 'Voltage Pulse Exact Fit', 'Current Pulse Exact Fit', 'Fontsize', 20, 'Location', 'eastoutside')
legend boxoff 
xlim([1,10])
f.Position = [100 100 910 400];
%HPPC EXAMPLE
export_fig fig25.png -transparent -m3 -r300



%now we plot the full cell with all deg mechanism


figure(6)
%copy paste reaction rate and Wfunc
semilogy(linspace(0,0.25,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(1:9,1),'LineWidth',2, 'Color', map_purple_less(end,:))
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(1:9,1),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(1:9,1),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
xlim([1,9])
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
export_fig fig26.png -transparent -m3 -r300


figure(7)
%copy paste reaction rate and Wfunc
plot(linspace(1,0.9,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,2),'LineWidth',2, 'Color', map_purple_less(end-4,:))
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(1:9,2),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,2),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
xlim([1,9])
export_fig fig27.png -transparent -m3 -r300


figure(8)
%copy paste reaction rate and Wfunc
plot(linspace(1000,900,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--', 'HandleVisibility','off')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'HandleVisibility','off')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(1:9,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'LineStyle', ':', 'HandleVisibility','off')
load('../current_pulse_results/full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,3)*1000,'LineWidth',2, 'Color', map_purple_less(end-7,:), 'LineStyle', '-.', 'HandleVisibility','off')
% load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_charge_linear.mat')
% plot(zPulse(:,3)*1000,'o','LineWidth',2, 'Color', map_purple_less(end-4,:))
% hold on
% load('../full_cell_pulse_fit_c_lyte_rescaled_together_full_cell_L2_discharge_linear.mat')
% plot(zPulse(:,3)*1000,'x','LineWidth',2, 'Color', map_purple_less(end-4,:))
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
xlim([1,9])
% legend('Discharge Fitted Solution' , 'Charge Fitted Solution', 'Fontsize', 20, 'Location', 'southwest')
% legend boxoff
%HPPC EXAMPLE
export_fig fig28.png -transparent -m3 -r300


figure(9)
%copy paste reaction rate and Wfunc
semilogy(linspace(0,0.25,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(1:9,4),'LineWidth',2, 'Color', map_purple_less(end,:))
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_exact.mat')
semilogy(zPulse(1:9,4),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
semilogy(zPulse(1:9,4),'LineWidth',2, 'Color', map_purple_less(end,:), 'LineStyle', '-.')

set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
% legend('Exact Solution', 'Linear Fitted Solution', 'Exact Fitted Solution')
%HPPC EXAMPLE
xlim([1,9])
export_fig fig29.png -transparent -m3 -r300

f = figure(10)
%copy paste reaction rate and Wfunc
plot(linspace(1,0.9,10),'LineWidth',2, 'Color', 'black', 'LineStyle', '--')
hold on
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:))
load('../full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_exact.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', ':')
load('../current_pulse_results/full_cell_pulse_fit_all_deg_rescaled_together_full_cell_L2_linear.mat')
plot(zPulse(1:9,5),'LineWidth',2, 'Color', map_purple_less(end-4,:), 'LineStyle', '-.')
set(gca,'fontsize',12,'LineWidth',2,'FontName','Helvetica')
%legend('Exact Solution', 'Voltage Pulse Linear Fit', 'Voltage Pulse Exact Fit', 'Current Pulse Exact Fit', 'Fontsize', 20, 'Location', 'southwest')
legend('Exact Solution', 'Voltage Pulse Linear Fit', 'Voltage Pulse Exact Fit', 'Current Pulse Exact Fit', 'Fontsize', 20, 'Location', 'eastoutside')
legend boxoff 
f.Position = [100 100 910 400];
xlim([1,9])
%HPPC EXAMPLE
export_fig fig30.png -transparent -m3 -r300



function [x_red, y_red] = get_red_points(eta, c, k0, lambda)
x_red = [];
y_red = [];
for i = 1:length(c)
    x_red = [x_red, -eta, eta];
    y_red = [y_red, dideta(c(i), -eta, k0, lambda), dideta(c(i), eta, k0, lambda)];
end

end

function [x_red, y_red] = get_red_points2(eta, c, lambda)
x_red = [];
y_red = [];
for i = 1:length(c)
    x_red = [x_red, -eta, eta];
    y_red = [y_red, iredoi(c(i), -eta, lambda), iredoi(c(i), eta, lambda)];
end

end



function out = dhelper_fundetaf(eta_f, lmbda)
%derivative of helper function for CIET
out = (eta_f.*exp(-(lmbda - (eta_f.^2 + lmbda.^(1/2) + 1).^(1/2)).^2./(4*lmbda))) ...
    ./((exp(-eta_f) + 1).*(eta_f.^2 + lmbda^(1/2) + 1).^(1/2)) - ...
    (lmbda^(1/2).*pi.^(1/2)*exp(-eta_f).*(erf((lmbda - (eta_f.^2 + ...
    lmbda^(1/2) + 1).^(1/2))./(2*lmbda^(1/2))) - 1))./(exp(-eta_f) + 1).^2;
end

function out = dideta(c, eta, k0, lambda)
%dideta function for CIET
% muh = mu_c(c, params);
% eta = muh - mures;
% etaf = eta - log(c);
etaf = eta - log(c);
out = k0.*(1-c)/sqrt(4*pi*lambda).*...
     (-dhelper_fundetaf(-etaf, lambda) - c.*dhelper_fundetaf(etaf, lambda));
end

function out = i(c, eta, k0, lambda)
%dideta function for CIET
% muh = mu_c(c, params);
% eta = muh - mures;
% etaf = eta - log(c);
etaf = eta - log(c);
out = k0.*(1-c)/sqrt(4*pi*lambda).*...
     (helper_fun(-etaf, lambda) - c.*helper_fun(etaf, lambda));
end



function rxn = iredoi(c, eta, lambda)
%ired/i fraction used in electrolyte fitness value R is dc/dt, from the
%population balance equation. can be BV or linear
% %set symmetry coefs
%muh = mu_c(c, params_nondim);
%for simplicity in calculation
etaf = eta - log(c);
rxn = helper_fun(-etaf, lambda)./(helper_fun(-etaf, lambda)-c.*helper_fun(etaf, lambda));
end



function vq = interp_colors(map_less, len)
vq = zeros(len, 3);
for i = 1:3
    xq = linspace(1, size(map_less,1), len);
    vq(:,i) = interp1(1:size(map_less,1), map_less(:,i), xq);
end
end


function out = helper_fun(eta_f, lmbda)
out = (sqrt(pi*lmbda)./(1+exp(-eta_f)).*...
    (1-erf((lmbda-sqrt(1+sqrt(lmbda)+eta_f.^2))./(2*sqrt(lmbda)))));
end



