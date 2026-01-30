clear;
clc;
close all;

%% Select drug-specific mean parameters needed for simulations
% The mean values are presented in Figure 3a-d

drug = "Amp";
% drug = "Cip";
% drug = "Kan";

if drug == "Amp"
    f_1 = 10^-2.708;
    k_1 = 8.374;
    k_2 = 0.816;
    init_pers_perc_10000x = 10^(-1.512 - 1); % Minus 1 to get 10,000x init %P
elseif drug == "Cip"
    f_1 = 10^-3.691;
    k_1 = 15.185;
    k_2 = 0.258;
    init_pers_perc_10000x = 10^(-2.815 - 1); % Minus 1 to get 10,000x init %P
elseif drug == "Kan"
    f_1 = 10^-6.052;
    k_1 = 22.896;
    k_2 = 0.629;
    init_pers_perc_10000x = 10^(-5.379 - 1); % Minus 1 to get 10,000x init %P
end

%% Set up the rest of the simulation information

% Assume N_0 is 10^8.25 CFU/mL (a typical value, I'm just going to normalize it away in percent survival anyways)
N_0 = 10^8.25;

% Calculate P_0 from the initial persister percent
P_0 = N_0 * (init_pers_perc_10000x/100);

% Calculate the initial number of susceptible cells
S_0 = N_0 - P_0;

% Organize initial conditions
init_conds = [S_0, P_0];

% Set simulation time
t_span = [0,8];

%% Simulation

[t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k_1,k_2,f_1), t_span, init_conds);

%% Calculate total cell amounts
N_of_t = sum(y_ode,2);
P_of_t = y_ode(:,2);

%% Plot simulation dynamics

% Figure size
fig1_pixels_width=500;
fig1_pixels_height=417;

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];
ax1 = gca;

hold on;

% Plot the ODE simulation results (N and P)
plot(t_ode,(100*N_of_t/N_of_t(1)),'Color','k','LineStyle','-','LineWidth',3);
plot(t_ode,(100*P_of_t/N_of_t(1)),'Color',[0.3, 0.7, 0.3],'LineStyle','-','LineWidth',3);

title(drug)
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,3]);
xticks(0:0.5:3);

if drug=="Amp"
    ylim([0,0.025]);
    yticks(0:0.005:0.025);
end

grid on;
set(ax1,'GridAlpha',0.4);

legend("Total Cells","Persister Cells",'Location','northeast','orientation','vertical');

ax1.FontSize = 28;

% Save plot
fig1_plot_name = strcat(drug,"_MeanParamTimeKillPersisterDynamics_Zoomed_10000x.png");
saveas(f1,fig1_plot_name);


%% Model ODEs
function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1)-k_2*y(2); % dydt(2,1) is dP_di/dt
end