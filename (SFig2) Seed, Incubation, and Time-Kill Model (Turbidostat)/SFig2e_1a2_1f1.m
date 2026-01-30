% This code will simulate the following scenario using the following models

% Phase 1: Exponential Phase in Turbidostat - Asssume that all Type I 
% (stress-induced) persisters are diluted out of the turbidostat and the 
% fraction of Type II (spontaneous) persisters reaches an equilibrium.
    % Model persister population dynamcis via the combined WT model from 
    % Balaban 2004.

% Action 1: Dilution of exponential phase cells into fresh media at a range
% of dilution factors.

% Phase 2: Lag Phase & Exponential Phase
    % Model persister population dynamcis via the combined WT model from 
    % Balaban 2004.

% Action 2: Addition of antibiotics to the exponential phase culture.
% Assume that persister types can be combined (i.e., have the same 
% interactions and same parameters) during antibiotic treatment

% Phase 3: Treatment Phase (Time-Kill Assay)
    % Use our model for antibiotic-induced persister formation to create 
    % the time-kill assay curves.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clears
clear;
clc;

%% Exponential Phase in Turbidostat (Balaban 2004 WT model)

% Balaban WT model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_2 = 1.2*10^-6*(1/60); % Per cell per hour times (1h / 60min), Given in Balaban Fig 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b_1 = 0.05*(1/60); % Per cell per hour times (1h / 60min), Given in Balaban Fig 2
b_2 = 0.1*(1/60); % Per cell per hour times (1h / 60min), Given in Balaban Fig 2
% mu_n = 2.01*(1/60); % Per hour times (1h / 60min), Estimated from Balaban 2004 Figure 3B
mu_n = 1.53*(1/60); % Per hour times (1h / 60min), Calculated from our own data
mu_p1 = 0; % Per hour times (1h / 60min), Approximately zero
mu_p2 = 0.13*(1/60); % Per hour times (1h / 60min), Fitted to match Fig 3b in Balaban 2004

% Culture Parameters
exp_cells_tot = 2*10^8; % Total CFU/ml in turbidostat at equilibrium

% Simulate the turbidostat's approach to equilibrium
t_turb_sim = [0,12]*60; % h times (60min / 1h)
n_0_frac = 0.9995; % S_0
p_1_0_frac = 0.0005; % P_1_0, assumed
p_2_0_frac = 0; % P_2_0
init_conds = [n_0_frac, p_1_0_frac, p_2_0_frac]*exp_cells_tot;
[t_ode_turb, y_ode_turb] = ode15s(@(t,y) balaban_2004_ODEs_turbidostat(t,y,a_2,b_1,b_2,mu_n,mu_p1,mu_p2), t_turb_sim, init_conds);

% Note: Comment in this plotting code to see the turbidostat equilibrate over time.
% figure(11);
% ax11 = gca;
% semilogy(t_ode_turb/60, y_ode_turb(:,1),'r','LineWidth',2);
% hold on;
% semilogy(t_ode_turb/60, y_ode_turb(:,2),'c','LineWidth',2);
% semilogy(t_ode_turb/60, y_ode_turb(:,3),'b','LineWidth',2);
% title("Turbidostat Equilibration");
% xlabel("Time (h)");
% ylabel("CFU/mL");
% ylim([10^0,10^9]);
% grid on;
% legend("Susceptible Cells","Type I Persisters","Type II Persisters","Location","East");
% set(ax11,"FontSize",14);


% Results
S_turb = y_ode_turb(end,1); % Susceptible cells in exponential phase culture
P_1_turb = 0; % Assume that there are no Type I persister cells in the exponential phase culture
P_2_turb = y_ode_turb(end,3); % Type II persister cells in exponential phase culture
P_tot_turb = P_1_turb + P_2_turb;
N_tot_turb = S_turb + P_tot_turb;

%% Dilution into fresh media

% Model parameters
dil_facts = [20; 200; 2000]; % dilution factors

% Initialize storage arrays for simulation results
n_dilution_factors = length(dil_facts);
results_storage = cell(n_dilution_factors,1);
% Columns in each cell: time (minutes) | S | P_1 | P_2 | P_tot | N_tot

% Results
for i = 1:n_dilution_factors
    results_storage{i} = [0, S_turb/dil_facts(i), P_1_turb/dil_facts(i), P_2_turb/dil_facts(i), P_tot_turb/dil_facts(i), N_tot_turb/dil_facts(i)];
end

%% Lag phase and exponential phase (Balaban 2004 WT model)

%% Lag phase 

% Assume no growth for any cell type

t_lag = 5; % minutes

% Loop through simulations
if t_lag>0
    for i = 1:n_dilution_factors
        t_span = [0, t_lag]; % Time to simulate across (min)
        y_0 = results_storage{i}(end,2:4); % ODE solver initial conditions (for S, P1, P2)
        [t_temp, y_temp] = ode45(@(t,y) balaban_2004_ODEs(t,y,a_2,b_1,b_2,0,0,0), t_span, y_0);

        % Remove initial conditions from results (because they're already
        % saved)
        t_temp(1) = [];
        y_temp(1,:) = [];
    
        % Append results to storage
        new_tpts_temp = size(t_temp,1);
        results_storage{i}(end+1:end+new_tpts_temp,:) = [t_temp, y_temp, y_temp(:,2)+y_temp(:,3),sum(y_temp,2)];
    
    end
end


%% Exponential Phase

% For each culture, give the time between dilution and the start of the
% time-kill assay.
% Minus 5 min adjustment to better align N_0 for the data and model
t_treatment = [120, 210, 300]-5; % 100x, 1000x, 10000x

% Loop through simulations
for i = 1:n_dilution_factors
    t_span = [t_lag, t_treatment(i)]; % Time to simulate across (min)
    y_0 = results_storage{i}(end,2:4); % ODE solver initial conditions (for S, P1, P2)
    [t_temp, y_temp] = ode45(@(t,y) balaban_2004_ODEs(t,y,a_2,b_1,b_2,mu_n,mu_p1,mu_p2), t_span, y_0);

    % Remove initial conditions from results (because they're already
    % saved)
    t_temp(1) = [];
    y_temp(1,:) = [];

    % Append results to storage
    new_tpts_temp = size(t_temp,1);
    results_storage{i}(end+1:end+new_tpts_temp,:) = [t_temp, y_temp, y_temp(:,2)+y_temp(:,3),sum(y_temp,2)];

end

%% Addition of antibiotics

%% Treatment Phase (Time-Kill Assay, Our Drug-Induced Persistence Model)

% Assume that P_1 and P_2 collapse into a single P_tot population

% Our model parameters (from SDTK fitting of Amp_WT_Turb data)
k_1 = 6.757*(1/60); % deaths per cell per hour times (1h / 60min)
k_2 = 0.793*(1/60); % deaths per cell per hour times (1h / 60min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_1 = 2.464*(10^-3)*(1/60); % drug induction of persistence per cell per hour times (1h / 60min)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time spent in treatment phase
treat_time = 8*(60/1); % hours times (60min / 1h)

% Loop through simulations
for i = 1:n_dilution_factors
    t_span = [t_treatment(i), t_treatment(i)+treat_time]; % Time to simulate across (min)
    y_0 = results_storage{i}(end,[2,5]); % ODE solver initial conditions (for S, P_tot)
    [t_temp, y_temp] = ode45(@(t,y) our_ODEs(t,y,k_1,k_2,f_1), t_span, y_0);

    % Remove initial conditions from results (because they're already
    % saved)
    t_temp(1) = [];
    y_temp(1,:) = [];

    % Append results to storage
    new_tpts_temp = size(t_temp,1);
    results_storage{i}(end+1:end+new_tpts_temp,:) = [t_temp, y_temp(:,1), NaN*ones(new_tpts_temp,2), y_temp(:,2), y_temp(:,1)+y_temp(:,2)];

end

%% Load our experimental results for comparison

expdata = readmatrix("Amp_WT_Turb_042525_Data.txt");

t_expdata_plot_h = expdata(:,1)/60;

colors_RGB = [1.0, 0.0, 0.0; % red
              0.3, 0.7, 0.3; % green
              0.0, 0.0, 1.0]; % blue



%% Plot results with aligned antibiotic start time

% Size of figure window
scale = 200;
pixels_width=round(2.125*scale);
pixels_height=round(1.75*scale);

% Create figure window
f1=figure(1);
f1.Position = [0, 0, pixels_width, pixels_height];

hold on;

% Loop through plotting
for i = 1:n_dilution_factors
    t_plot_h = results_storage{i}(:,1)/60 - t_treatment(i)/60;
    log10_S_plot = log10(results_storage{i}(:,2));
    log10_P_1_plot = log10(results_storage{i}(:,3));
    log10_P_2_plot = log10(results_storage{i}(:,4));
    log10_P_tot_plot = log10(results_storage{i}(:,5));
    log10_N_tot_plot = log10(results_storage{i}(:,6));

    semilogy(t_plot_h,log10_S_plot,'b',"LineWidth",2);
    semilogy(t_plot_h,log10_P_1_plot,'Color',[12,122,6]/255,"LineWidth",6);
    semilogy(t_plot_h,log10_P_2_plot,'Color',[141,8,194]/255,"LineWidth",6);
    semilogy(t_plot_h,log10_P_tot_plot,'r',"LineWidth",2);
    semilogy(t_plot_h,log10_N_tot_plot,'k',"LineWidth",2);

    scatter(t_expdata_plot_h,expdata(:,1+i),25,'MarkerEdgeColor',colors_RGB(i,:),'MarkerFaceColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);
    errorbar(t_expdata_plot_h,expdata(:,1+i),expdata(:,1+n_dilution_factors+i),'o','color','k');
end

xlabel('Time (h)');
ylabel('log_{10}(CFU/mL)');

legend('S','','P_2','P_{tot}','N');

title('a_2 = Reported Value   &   f_1 = This Study');

grid on;

ylim([1,8.5]);
yticks(1:1:8);
xlim([-5,8]);
xticks(-5:1:8);

% Font Size
axis1 = gca; % Get the axis object for the plots
axis1.FontSize=18; % Change the font size


% Save figure
saveas(f1,"SFig2e_1a2_1f1.png");


%% ODE functions

%% balaban_2004_ODEs
function dydt = balaban_2004_ODEs_turbidostat(~,y,a_2,b_1,b_2,mu_n,mu_p1,mu_p2)
    
    % Extract states from y
    n = y(1); % Equivalent to S
    p_1 = y(2);
    p_2 = y(3);

    % Calculate dilution rate (weighted average of growth rates)
    D = (n*mu_n + p_1*mu_p1 + p_2*mu_p2)/(n+p_1+p_2);

    % Initialize dydt
    dydt = zeros(length(y),1);

    % Calculate derivatives
    dydt(1) = -1*a_2*n + b_1*p_1 + b_2*p_2 + mu_n*n - D*n; % dn/dt
    dydt(2) = -1*b_1*p_1 + mu_p1*p_1 - D*p_1; % dp1/dt
    dydt(3) = a_2*n - b_2*p_2 + mu_p2*p_2 - D*p_2; % dp2/dt

end

%% balaban_2004_ODEs
function dydt = balaban_2004_ODEs(~,y,a_2,b_1,b_2,mu_n,mu_p1,mu_p2)
    
    % Extract states from y
    n = y(1); % Equivalent to S
    p_1 = y(2);
    p_2 = y(3);

    % Initialize dydt
    dydt = zeros(length(y),1);

    % Calculate derivatives
    dydt(1) = -1*a_2*n + b_1*p_1 + b_2*p_2 + mu_n*n; % dn/dt
    dydt(2) = -1*b_1*p_1 + mu_p1*p_1; % dp1/dt
    dydt(3) = a_2*n - b_2*p_2 + mu_p2*p_2; % dp2/dt

end

%% our_ODEs
function dydt = our_ODEs(~,y,k_1,k_2,f_1)
    
    % Extract states from y
    S = y(1);
    P = y(2);

    % Initialize dydt
    dydt = zeros(length(y),1);

    % Calculate derivatives
    dydt(1) = -1*k_1*S - f_1*S; % dS/dt
    dydt(2) = f_1*S - k_2*P; % dP/dt

end