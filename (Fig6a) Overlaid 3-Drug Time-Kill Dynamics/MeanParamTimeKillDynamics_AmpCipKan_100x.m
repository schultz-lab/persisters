%% Clear
clear;
clc;
close all;

%% Collect the drug-specific mean parameters needed for simulations
% The mean values are presented in Figure 3a-d

% Amp mean parameters
f_1_Amp = 10^-2.708;
k_1_Amp = 8.374;
k_2_Amp = 0.816;
init_pers_perc_1000x_Amp = 10^-1.512;

% Cip mean parameters
f_1_Cip = 10^-3.691;
k_1_Cip = 15.185;
k_2_Cip = 0.258;
init_pers_perc_1000x_Cip = 10^-2.815;

% Kan mean parameters
f_1_Kan = 10^-6.052;
k_1_Kan = 22.896;
k_2_Kan = 0.629;
init_pers_perc_1000x_Kan = 10^-5.379;


%% Set the initial amounts of each population

% Assume N_0 is 10^8.25 CFU/mL (a typical value, I'm just going to normalize it away in percent survival anyways)
N_0 = 10^8.25;

% Convert to initial persister fraction for 100x
ip_frac_Amp100xdilute = 10*init_pers_perc_1000x_Amp/100;
ip_frac_Cip100xdilute = 10*init_pers_perc_1000x_Cip/100;
ip_frac_Kan100xdilute = 10*init_pers_perc_1000x_Kan/100;

% Calculate initial conditions
P_0_Amp = N_0*ip_frac_Amp100xdilute;
S_0_Amp = N_0 - P_0_Amp;
P_0_Cip = N_0*ip_frac_Cip100xdilute;
S_0_Cip = N_0 - P_0_Cip;
P_0_Kan = N_0*ip_frac_Kan100xdilute;
S_0_Kan = N_0 - P_0_Kan;


%% Simulate the dynamics

% Simulation time
t_span = [0,8];

% Amp treatment
init_conds = [S_0_Amp, P_0_Amp];
[t_ode_Amp, y_ode_Amp] = ode45(@(t,y) model_ODEs(t,y,k_1_Amp,k_2_Amp,f_1_Amp), t_span, init_conds);

% Cip treatment
init_conds = [S_0_Cip, P_0_Cip];
[t_ode_Cip, y_ode_Cip] = ode45(@(t,y) model_ODEs(t,y,k_1_Cip,k_2_Cip,f_1_Cip), t_span, init_conds);

% Kan treatment
init_conds = [S_0_Kan, P_0_Kan];
[t_ode_Kan, y_ode_Kan] = ode45(@(t,y) model_ODEs(t,y,k_1_Kan,k_2_Kan,f_1_Kan), t_span, init_conds);

% Calculate total viable cells and take the log10 of the results
log10_SPN_Amp = log10([y_ode_Amp(:,1), y_ode_Amp(:,2), sum(y_ode_Amp,2)]);
log10_SPN_Cip = log10([y_ode_Cip(:,1), y_ode_Cip(:,2), sum(y_ode_Cip,2)]);
log10_SPN_Kan = log10([y_ode_Kan(:,1), y_ode_Kan(:,2), sum(y_ode_Kan,2)]);


%% Plot results

% Figure size
fig1_pixels_width=600*(5/6);
fig1_pixels_height=500*(5/6);

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];
ax1 = gca;

hold on;
% Plot the Amp treatment results
% N
plot(t_ode_Amp,log10_SPN_Amp(:,3)-log10_SPN_Amp(1,3)+2,'Color',[76, 176, 9]/255,'LineStyle','-','LineWidth',3);
% P_Amp
plot(t_ode_Amp,log10_SPN_Amp(:,2)-log10_SPN_Amp(1,3)+2,'Color',[76, 176, 9]/255,'LineStyle','--','LineWidth',3);

% Plot the Cip treatment results
% N
plot(t_ode_Cip,log10_SPN_Cip(:,3)-log10_SPN_Cip(1,3)+2,'Color',[7, 60, 240]/255,'LineStyle','-','LineWidth',3);
% P_Kan
plot(t_ode_Cip,log10_SPN_Cip(:,2)-log10_SPN_Cip(1,3)+2,'Color',[7, 60, 240]/255,'LineStyle','--','LineWidth',3);

% Plot the Kan treatment results
% N
plot(t_ode_Kan,log10_SPN_Kan(:,3)-log10_SPN_Kan(1,3)+2,'Color',[224, 111, 7]/255,'LineStyle','-','LineWidth',3);
% P_Kan
plot(t_ode_Kan,log10_SPN_Kan(:,2)-log10_SPN_Kan(1,3)+2,'Color',[224, 111, 7]/255,'LineStyle','--','LineWidth',3);

xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
ylim([-8,2.2]);
yticks(-8:1:2);
yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);

grid on;
set(ax1,'GridAlpha',0.4);
legend("Amp","","Cip","","Kan",'Location','northeast','orientation','vertical');

ax1.FontSize = 28;

% Save plot
fig1_plot_name = strcat("MeanParamTimeKillDynamics_AmpCipKan_100x.png");
saveas(f1,fig1_plot_name);


%% Functions

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end