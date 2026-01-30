%% Notes
% Whenever S (susceptible cells) and P (persister cells) are together in
% the same matrix, we will order them so S is #1 and P is #2.

%% Housekeeping
clear;
clc;
close all;

%% Set up ODE solver

t_0 = 0; % Start time of the simulation, unit: hours
t_end = 7; % End time of the simulation, unit: hours
t_span = [t_0, t_end];

N_0 = 10^8.25; % Assumed, a typical value that I'm just going to normalize
               % away in percent survival anyways)
init_pers_frac = 10^-5; % Assumed, a moderate value for 1000x culture
P_0 = N_0*init_pers_frac; % Initial persister cell population
S_0 = N_0 - P_0; % Initial susceptible cell population
init_conds = [S_0, P_0];

%%% Co-Looped Parameters %%%

%%%%%%% Rounded values from fits to FigS4b data %%%%%%%
k_1 = [11, 18, 29]; % Death rate of susceptible cells, unit: 1/h
k_2 = [1.3,0.7,0.9]; % Death rate of persister cells, unit: 1/h

n_cults = length(k_1);
%%% Co-Looped Parameters %%%

% No drug-induced persistence
f_1 = 0; % Rate of susceptible to persister induction, unit: 1/h.

%% Prepare storage of results
t_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
log10_SPN_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
                                       % each array entry is 3 columns: S, P, and N.

%% Loop ODE Solver Runs & Results Storage

for i = 1:n_cults

    % Run simulation
    [t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k_1(i),k_2(i),f_1), t_span, init_conds);

    % Store data
    t_ODEsolver{i} = t_ode;
    log10_SPN_ODEsolver{i} = log10([y_ode(:,1), y_ode(:,2), y_ode(:,1)+y_ode(:,2)]);

end

%% Plot Results

% Size of figure window
pixels_width=600;
pixels_height=500;

% Create figure window
f1=figure(1);
f1.Position = [0, 0, pixels_width, pixels_height];

% Define a color palette for plots
colors_RGB = [ 76, 176,   9;     % Amp
                7,  60, 240;     % Cip
              224, 111,   7]/255;     % Kan

hold on;

% Plot the ODE simulation results (N)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,3)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle','-','LineWidth',3);
end

% Plot the ODE simulation results (P)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,2)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',3);
end

% Figure Labels
text(0.05, 2.05, "N_0","Color",'k','FontSize',22,'FontWeight','bold');
text(0.05, -3.5, "P_0","Color",'k','FontSize',22,'FontWeight','bold');

% Generate legend entries
ax_legend = cell(n_cults,1);
ax_legend{1} = "1,000x Seed Dilution (Amp)";
ax_legend{2} = "1,000x Seed Dilution (Cip)";
ax_legend{3} = "1,000x Seed Dilution (Kan)";

legend(ax_legend,'Location','northeast','orientation','vertical');


xlabel("Time (h)");
ylabel("Percent Survival (%)");
grid on;
set(get(f1,'CurrentAxes'),'GridAlpha',0.4);
ylim([-7,2.5]);
yticks(-7:1:2);
yticklabels(["10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
xlim([0,7]);
xticks(0:1:7);

% Font Size
axis1 = gca;
axis1.FontSize=28;

saveas(f1,"SameCulure_DifDrugs_f1_equals_0.png");

%% model_ODEs function

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P 
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end
