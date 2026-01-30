clear;
clc;
close all;

%% Note: Comment in the drug for which you want to do the relative persister fraction analysis
% drug = "Amp";
% drug = "Cip";
% drug = "Kan";

%% Note: Set the timepoint at which you want to do the relative persister fraction analysis
t_frac = 3; % hours

%% Load SDTK fitting results
results_table_holder = load('FittingResultsSuperTable.mat');
results_table = results_table_holder.supertable;

%% Select experimental datasets

if drug == "Amp"
    experiment_indices = [1, 2, 3];
elseif drug == "Cip"
    experiment_indices = [4, 5, 6, 7];
elseif drug == "Kan"
    experiment_indices = [8, 9, 10, 11];
end
                                                            % Experiment
                                                            % Indices
datafile_list = ["SDTK_Amp_WT_Rep1_022625_Data.txt";        % 1  = Amp_WT_Rep1
                 "SDTK_Amp_WT_Rep2_040125_Data.txt";        % 2  = Amp_WT_Rep2
                 "SDTK_Amp_WT_Rep3_061025_Data.txt";        % 3  = Amp_WT_Rep3
                 "SDTK_Cip_WT_Rep1_021225_Data.txt";        % 4  = Cip_WT_Rep1
                 "SDTK_Cip_WT_Rep2_092625_Data.txt";        % 5  = Cip_WT_Rep2
                 "SDTK_Cip_WT_Rep3_100225_Data.txt";        % 6  = Cip_WT_Rep3
                 "SDTK_Cip_WT_Rep4_110525_Data.txt";        % 7  = Cip_WT_Rep4
                 "SDTK_Kan_WT_Rep1_051325_Data.txt";        % 8  = Kan_WT_Rep1
                 "SDTK_Kan_WT_Rep2_053025_Data.txt";        % 9  = Kan_WT_Rep2
                 "SDTK_Kan_WT_Rep3_100525_Data.txt";        % 10 = Kan_WT_Rep3
                 "SDTK_Kan_WT_Rep4_102625_Data.txt";        % 11 = Kan_WT_Rep4
                 "SDTK_Amp_WT_Turb_042525_Data.txt";        % 12 = Amp_WT_Turb
                 "SDTK_Amp_ppk_Rep1_060925_Data.txt";       % 13 = Amp_ppk_Rep1
                 "SDTK_Amp_ppk_Rep2_061225_Data.txt";       % 14 = Amp_ppk_Rep2
                 "SDTK_Amp_ppk_Rep3_102825_Data.txt";       % 15 = Amp_ppk_Rep3
                 "SDTK_Amp_relAspoT_Rep1_081925_Data.txt";  % 16 = Amp_relAspoT_Rep1
                 "SDTK_Amp_relAspoT_Rep2_090525_Data.txt";  % 17 = Amp_relAspoT_Rep2
                 "SDTK_Amp_relAspoT_Rep3_091125_Data.txt";  % 18 = Amp_relAspoT_Rep3
                 "SDTK_Amp_relAspoT_Rep4_101625_Data.txt"]; % 19 = Amp_relAspoT_Rep4

%% Loop through dataset analyses
n_experiments = length(experiment_indices);

Ppe_frac = zeros(n_experiments,3); % Pre-existing persister fraction, 3 columns for 100x, 1000x, 10000x

for m = 1:n_experiments
    experiment_index = experiment_indices(m);

    %% Load & extract experimental data
    
    % Load data
    datafile = datafile_list(experiment_index);
    datafile_path = strcat("Data_and_Descriptions/",datafile);
    data = readmatrix(datafile_path);
    
    % Convert time to hours
    data(:,1) = data(:,1)/60;
    
    % Determine number of cultures
    n_cults = (size(data,2)-1)/2; % Look at number of columns. 
                                  % Minus 1 for the time, /2 for mean & stddev
    
    % Determine index of last timepoint in each culture
    last_time_index = -1*ones(n_cults,1);
    for i = 1:n_cults
        last_time_index(i) = find(~isnan(data(:,1+i)),1,'last');
    end
    
    % Create a cell array of (time, mean cell ct, std dev cell ct) datasets for
    % each culture (100x, 1000x, 10000x)
    t_mean_stddev = {};
    for i = 1:n_cults
        t_mean_stddev{i} = data(1:last_time_index(i),[1,1+i,1+n_cults+i]);
    end

    %% Compile information needed for simulations

    N_0_vect = 10.^[t_mean_stddev{1}(1,2), t_mean_stddev{2}(1,2), t_mean_stddev{3}(1,2)]; % 100x, 1000x, 10000x

    init_pers_frac_vect = results_table.init_pers_frac_leastdilute_fitted(experiment_index)./[1, 10, 100];

    k_1 = results_table.k1_est_mean(experiment_index);

    k_2 = results_table.k2_est_mean(experiment_index);

    f_1 = results_table.f1_fitted(experiment_index);

    t_span = [0,8];

    %% Simulation
    
    % Create storage arrays for ODE solver solutions
    ODE_soln = cell(n_cults,1); % Save the general form of the ODE soln so I can use deval() later
    
    % Loop through ODE solver calls
    for i = 1:n_cults
        % Calculate initial conditions
        N_0 = N_0_vect(i); % Total cells
        P_pe_0 = N_0*init_pers_frac_vect(i); % Pre-existing persister cells
        P_di_0 = 0; % Drug-induced persisters
        S_0 = N_0-P_pe_0; % Susceptible cells
    
        % Simulation
        init_conds = [S_0, P_pe_0, P_di_0];
        ODE_soln{i} = ode45(@(t,y) model_ODEs_separateP(t,y,k_1,k_2,f_1), t_span, init_conds);
    end
    
    %% Calculate the persister fractions at a given time (t_frac)
    
    y_at_t_frac = cell(n_cults,1); % To store the ODE solution evaluated at 
                                   % t_frac for the 3 dilution conditions.
    for i = 1:n_cults
        y_at_t_frac{i} = deval(ODE_soln{i},t_frac);
    end
    
    % Calculate persister fractions
    for i = 1:n_cults
        Ppe_frac(m,i) = y_at_t_frac{i}(2)/(y_at_t_frac{i}(2)+y_at_t_frac{i}(3));
    end

end


%% Calculate the mean and std deviation of persister fractions
Ppe_frac_mean = mean(Ppe_frac,1);
Ppe_frac_stddev = std(Ppe_frac,0,1);


%% Make a stacked bar plot of the persister percents at t_frac

% Figure size
fig1_pixels_width=600;
fig1_pixels_height=400;

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];

hold on;

% Format results
bar_data = 100*[1-Ppe_frac_mean; Ppe_frac_mean]; % Units: %

% Plot results
ax1 = gca;

bar_labels = ["100x", "1,000x", "10,000x"];
bar_handle = bar(bar_labels,bar_data,'stacked','LineWidth',1);
set(bar_handle,'FaceColor','flat');
bar_handle(1).CData = [0,0,0];
bar_handle(2).CData = [0.8,0.8,0.8];

er=errorbar(categorical(bar_labels),100-bar_data(2,:),100*Ppe_frac_stddev,"LineWidth",5,"CapSize",25); % Units: %
er.Color = [0.8 0 0];                            
er.LineStyle = 'none';

xlabel("Seed Dilution");
ylabel("% of Total Persisters");
ylim([0,105]);
yticks(0:25:100);
h=legend("Drug-Induced","Pre-Existing","Orientation","horizontal","Location","northoutside");

ax1.YGrid = 'on';
set(ax1,'GridAlpha',0.5);
set(ax1,'GridLineWidth',1);
ax1.FontSize = 32;

% Save plot
saveas(f1,sprintf("Mean_%s_relative_persister_percents_at_%dh.png",drug,t_frac));


%% Model ODEs
function dydt = model_ODEs_separateP(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P_pe (pre-existing persisters)
    % y(3) is P_di (drug-induced persisters)
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = -k_2*y(2); % dydt(2,1) is dP_pe/dt
    dydt(3,1) = f_1*y(1)-k_2*y(3); % dydt(3,1) is dP_di/dt
end
