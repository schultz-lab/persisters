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
    init_pers_perc_10000x = 10^(-1.512 - 1); % Minus 1 to get 10,000x %P
elseif drug == "Cip"
    f_1 = 10^-3.691;
    k_1 = 15.185;
    k_2 = 0.258;
    init_pers_perc_10000x = 10^(-2.815 - 1); % Minus 1 to get 10,000x %P
elseif drug == "Kan"
    f_1 = 10^-6.052;
    k_1 = 22.896;
    k_2 = 0.629;
    init_pers_perc_10000x = 10^(-5.379 - 1); % Minus 1 to get 10,000x %P
end

%% Set up the rest of the simulation information

% Assume N_0 is 10^8.25 CFU/mL (a typical value, I'm just going to normalize it away in percent survival anyways)
N_0 = 10^8.25;

% Calculate P_pe_0 from the initial persister percent
P_pe_0 = N_0 * (init_pers_perc_10000x/100);

% The initial number of drug-induced persisters is 0
P_di_0 = 0;

% Calculate the initial number of susceptible cells
S_0 = N_0 - P_pe_0 - P_di_0;

% Organize initial conditions
init_conds = [S_0, P_pe_0, P_di_0];

% Set simulation time
t_span = [0,8];

%% Simulation - For extracting data at specific timepoints (Fig S3c)
% NOTE: Comment in this simulation section, comment out all of the code
% after it (except the function), and select a t_eval timepoint to extract
% the data shown in Fig 5c.

% t_eval = 5*(1/60); % minutes times (1h / 60min)
% 
% ode_soln = ode45(@(t,y) model_ODEs_separateP(t,y,k_1,k_2,f_1), t_span, init_conds);
% 
% y_ode_at_t_eval = deval(ode_soln, t_eval);
% 
% P_tot_at_t_eval = sum(y_ode_at_t_eval(2:3));
% 
% percent_Ppe_at_t_eval = 100*y_ode_at_t_eval(2)/P_tot_at_t_eval;
% percent_Pdi_at_t_eval = 100*y_ode_at_t_eval(3)/P_tot_at_t_eval;
% 
% percents_of_interest_at_t_eval = [percent_Ppe_at_t_eval, percent_Pdi_at_t_eval];

%% Simulation - For making stacked plots of persister formation dynamics (Fig S3b)

[t_ode,y_ode] = ode45(@(t,y) model_ODEs_separateP(t,y,k_1,k_2,f_1), t_span, init_conds);

%% Calculate Percents
N_of_t = sum(y_ode,2);

percent_S = 100*y_ode(:,1)./N_of_t;
percent_Ppe = 100*y_ode(:,2)./N_of_t;
percent_Pdi = 100*y_ode(:,3)./N_of_t;

%% Plot simulation dynamics as stacked percentages

% Figure size
fig1_pixels_width=600;
fig1_pixels_height=500;

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];
ax1 = gca;

% Format data for area plot (each row is a timepoint)
area_plot_data = [percent_Pdi, percent_Ppe, percent_S];

% Area plot
area(t_ode,area_plot_data);

xlabel("Time (h)");
ylabel("Percent of Viable Cells");
xlim([0,8]);
xticks(0:8);
ylim([0,100]);
yticks(0:25:100);
legend('Drug-Induced Persisters','Pre-existing Persisters','Susceptible Cells','Location','southeast','orientation','vertical');
grid on;
set(ax1,'GridAlpha',0.6);
ax1.GridLineWidth = 1;
set(ax1,'Layer','top');
ax1.FontSize = 28;

% Save plot
fig1_plot_name = strcat(drug,"_MeanParamTimeKillPersisterDynamics_Stacked_Amp.png");
saveas(f1,fig1_plot_name);


%% Model ODEs
function dydt = model_ODEs_separateP(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P_pe (pre-existing persisters)
    % y(3) is P_di (drug-induced persisters)
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = -k_2*y(2); % dydt(2,1) is dP_pe/dt
    dydt(3,1) = f_1*y(1)-k_2*y(3); % dydt(3,1) is dP_di/dt
end
