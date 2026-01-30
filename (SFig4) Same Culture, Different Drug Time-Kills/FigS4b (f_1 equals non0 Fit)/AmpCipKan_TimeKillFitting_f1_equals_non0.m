%% Notes


%% Clear
clear;
clc;
close all; % figures

%% Select Dataset to Analyze

datafile = "SDTK_SameCult_DifDrugs_111925_Data.txt";

experiment_nickname = "SameCult_DifDrugs";


%% Load information about the selected dataset

% What the end time (h) of the k_1 linear region is.
t_k1_end = [0.5 0.5 0.5]; % Amp, Cip, Kan
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit

% What the start time (h) of the k_2 linear region is
t_k2_start = [2 1 1]; % Amp, Cip, Kan
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit

% Set of dilution factors used
dil_fact_set = [1000, 1000, 1000]; % Amp, Cip, Kan

% First Culture Type
first_cult_type = "batch";


%% Define a color palette for plots
colors_RGB = [ 76, 176,   9;     % Amp
                7,  60, 240;     % Cip
              224, 111,   7]/255;     % Kan


%% Load Dataset & Extract Data

% Load data
data = readmatrix(datafile);

% Convert time to hours
data(:,1) = data(:,1)/60;

% Determine number of cultures
n_cults = (size(data,2)-1)/2; % Look at number of columns. 
                              % Minus 1 for the time, /2 for mean & stddev

% Create a cell array of (time, mean cell ct, std dev cell ct) datasets for
% each dilution condition
t_mean_stddev = cell(n_cults,1);
t_mean_stddev{1} = data(:,[1,2,5]); % Amp 1000x
t_mean_stddev{2} = data(:,[1,3,6]); % Cip 1000x
t_mean_stddev{3} = data(:,[1,4,7]); % Kan 1000x


%% Plot Data (Raw)

% Figure size
fig1_pixels_width=600;
fig1_pixels_height=500;

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];
ax1 = gca;

hold on;
% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title("Time-Kill Assay Data");
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,7]);
xticks(0:7);
ylim([-0.3,9]);
yticks(0:9);
grid on;

% Generate legend entries
ax1_legends = cell(n_cults,1);
ax1_legends{1} = strcat("1,000x Seed Dilution (Amp)");
ax1_legends{2} = strcat("1,000x Seed Dilution (Cip)");
ax1_legends{3} = strcat("1,000x Seed Dilution (Kan)");

legend(ax1_legends,'Location','northeast','orientation','vertical');
ax1.FontSize = 14;

% Save plot
fig1_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot1_ExperimentalData.png");
saveas(f1,fig1_plot_name);


%% Plot Data (Normalized to 1)

% Figure size
fig2_pixels_width=600;
fig2_pixels_height=500;

% Create figure and axis
f2=figure(2);
f2.Position = [0, 0, fig2_pixels_width, fig2_pixels_height];
ax2 = gca;

hold on;
% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title("Time-Kill Assay Data");
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,7]);
xticks(0:7);
ylim([-7,2.2]);
yticks(-7:1:2);
yticklabels(["10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
grid on;

% Generate legend entries
ax2_legends = cell(n_cults,1);
ax2_legends{1} = strcat("1,000x Seed Dilution (Amp)");
ax2_legends{2} = strcat("1,000x Seed Dilution (Cip)");
ax2_legends{3} = strcat("1,000x Seed Dilution (Kan)");

legend(ax2_legends,'Location','northeast','orientation','vertical');
ax2.FontSize = 14;

% Save plot
fig2_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot1normed_ExperimentalData.png");
saveas(f2,fig2_plot_name);



%% Estimate k_1 from the first linear region

% Truncate data to the first linear region
t_mean_stddev_k1est = cell(size(t_mean_stddev));
for i = 1:n_cults
    n_tpts = sum(t_mean_stddev{i}(:,1) <= t_k1_end(i));
    t_mean_stddev_k1est{i} = t_mean_stddev{i}(1:n_tpts,:);
end

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k1_fits = zeros(3,n_cults);
% STRUCTURE OF INFO: for each column (culture), 
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_1 estimate
linear_k1fit_coefs_storage = zeros(n_cults,2);
% Structure: each row (culture) contains the 2 polynomial coefficients of
% the linear fit to k_1. These are stored for later plotting.

% Loop through linear fits and conversion to k_1 estimates
for i = 1:n_cults
    linear_fit_coefs = polyfit(t_mean_stddev_k1est{i}(:,1),t_mean_stddev_k1est{i}(:,2),poly_order);
    linear_k1fit_coefs_storage(i,:)=linear_fit_coefs;

    r_sqrd_mat = corrcoef(t_mean_stddev_k1est{i}(:,2), polyval(linear_fit_coefs,t_mean_stddev_k1est{i}(:,1)));
    k1_fits(1,i) = r_sqrd_mat(1,2); % R_squared

    E = t_mean_stddev_k1est{i}(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k1est{i}(:,1));
    SSE = sum(E.^2);
    k1_fits(2,i) = sqrt(SSE/length(t_mean_stddev_k1est{i}(:,2))); % RMSE

    % Multiply the slopes of the linear fits by ln(10) to estimate k_1 values
    k1_fits(3,i) = -1*linear_fit_coefs(1)*log(10); % estimated k_1
end

% Extract the estimated k_1 value
k1_amp_est = k1_fits(3,1);
k1_cip_est = k1_fits(3,2);
k1_kan_est = k1_fits(3,3);

k1_ests = [k1_amp_est; k1_cip_est; k1_kan_est];


%% Plot k_1 fit (Raw)

% Figure size
fig3_pixels_width=600;
fig3_pixels_height=500;

% Create figure and axis
f3=figure(3);
f3.Position = [0, 0, fig3_pixels_width, fig3_pixels_height];
ax3 = gca;

hold on;
% Plot the linear fits as dotted lines
for i = 1:n_cults
    plot(0:0.01:8,polyval(linear_k1fit_coefs_storage(i,:),0:0.01:8),'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title("Fitting for k_1 Estimation");
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,7]);
xticks(0:7);
ylim([-0.3,9]);
yticks(0:9);
grid on;

% Generate legend entries
ax3_legends = cell(n_cults,1);
ax3_legends{1} = strcat("1,000x Seed Dilution (Amp) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,1)));
ax3_legends{2} = strcat("1,000x Seed Dilution (Cip) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,2)));
ax3_legends{3} = strcat("1,000x Seed Dilution (Kan) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,3)));


legend(ax3_legends,'Location','northeast','orientation','vertical');
ax3.FontSize = 14;

% Save plot
fig3_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot2_k1Estimation.png");
saveas(f3,fig3_plot_name);


%% Plot k_1 fit (Normalized to 1)

% Figure size
fig4_pixels_width=600;
fig4_pixels_height=500;

% Create figure and axis
f4=figure(4);
f4.Position = [0, 0, fig4_pixels_width, fig4_pixels_height];
ax4 = gca;

hold on;
% Plot the linear fits as dotted lines
for i = 1:n_cults
    plot(0:0.01:8,polyval(linear_k1fit_coefs_storage(i,:),0:0.01:8)-t_mean_stddev{i}(1,2)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title("Fitting for k_1 Estimation");
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,7]);
xticks(0:7);
ylim([-7,2.2]);
yticks(-7:1:2);
yticklabels(["10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
grid on;

% Generate legend entries
ax4_legends = cell(n_cults,1);
ax4_legends{1} = strcat("1,000x Seed Dilution (Amp) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,1)));
ax4_legends{2} = strcat("1,000x Seed Dilution (Cip) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,2)));
ax4_legends{3} = strcat("1,000x Seed Dilution (Kan) k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,3)));

legend(ax4_legends,'Location','northeast','orientation','vertical');
ax4.FontSize = 14;

% Save plot
fig4_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot2normed_k1Estimation.png");
saveas(f4,fig4_plot_name);



%% Estimate k_2 from the second linear region

% Truncate data to the second linear region
t_mean_stddev_k2est = cell(size(t_mean_stddev));
for i = 1:n_cults
    n_tpts = sum(t_mean_stddev{i}(:,1) >= t_k2_start(i));
    t_mean_stddev_k2est{i} = t_mean_stddev{i}((end-n_tpts+1):end,:);
end

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k2_fits = zeros(3,n_cults);
% STRUCTURE OF INFO: for each column (culture), 
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_1 estimate
linear_k2fit_coefs_storage = zeros(n_cults,2);
% Structure: each row (culture) contains the 2 polynomial coefficients of
% the linear fit to k_2. These are stored for later plotting.

% Loop through linear fits and conversion to k_2 estimates
for i = 1:n_cults
    linear_fit_coefs = polyfit(t_mean_stddev_k2est{i}(:,1),t_mean_stddev_k2est{i}(:,2),poly_order);
    linear_k2fit_coefs_storage(i,:)=linear_fit_coefs;

    r_sqrd_mat = corrcoef(t_mean_stddev_k2est{i}(:,2), polyval(linear_fit_coefs,t_mean_stddev_k2est{i}(:,1)));
    k2_fits(1,i) = r_sqrd_mat(1,2); % R_squared

    E = t_mean_stddev_k2est{i}(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k2est{i}(:,1));
    SSE = sum(E.^2);
    k2_fits(2,i) = sqrt(SSE/length(t_mean_stddev_k2est{i}(:,2))); % RMSE

    % Multiply the slopes of the linear fits by ln(10) to estimate k_2 values
    k2_fits(3,i) = -1*linear_fit_coefs(1)*log(10); % estimated k_2
end

% Calculate the mean estimated k_2 value and its stddev
k2_amp_est = k2_fits(3,1);
k2_cip_est = k2_fits(3,2);
k2_kan_est = k2_fits(3,3);

k2_ests = [k2_amp_est; k2_cip_est; k2_kan_est];


%% Plot k_2 fit (Raw)

% Figure size
fig5_pixels_width=600;
fig5_pixels_height=500;

% Create figure and axis
f5=figure(5);
f5.Position = [0, 0, fig5_pixels_width, fig5_pixels_height];
ax5 = gca;

hold on;
% Plot the linear fits as dotted lines
for i = 1:n_cults
    plot(0:0.01:8,polyval(linear_k2fit_coefs_storage(i,:),0:0.01:8),'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title("Fitting for k_2 Estimation");
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,7]);
xticks(0:7);
ylim([-0.3,9]);
yticks(0:9);
grid on;

% Generate legend entries
ax5_legends = cell(n_cults,1);
ax5_legends{1} = strcat("1,000x Seed Dilution (Amp) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,1)));
ax5_legends{2} = strcat("1,000x Seed Dilution (Cip) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,2)));
ax5_legends{3} = strcat("1,000x Seed Dilution (Kan) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,3)));

legend(ax5_legends,'Location','northeast','orientation','vertical');
ax5.FontSize = 14;

% Save plot
fig5_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot3_k2Estimation.png");
saveas(f5,fig5_plot_name);


%% Plot k_2 fit (Normalized to 1)

% Figure size
fig6_pixels_width=600;
fig6_pixels_height=500;

% Create figure and axis
f6=figure(6);
f6.Position = [0, 0, fig6_pixels_width, fig6_pixels_height];
ax6 = gca;

hold on;
% Plot the linear fits as dotted lines
for i = 1:n_cults
    plot(0:0.01:8,polyval(linear_k2fit_coefs_storage(i,:),0:0.01:8)-t_mean_stddev{i}(1,2)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title("Fitting for k_2 Estimation");
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,7]);
xticks(0:7);
ylim([-7,2.2]);
yticks(-7:1:2);
yticklabels(["10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
grid on;

% Generate legend entries
ax6_legends = cell(n_cults,1);
ax6_legends{1} = strcat("1,000x Seed Dilution (Amp) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,1)));
ax6_legends{2} = strcat("1,000x Seed Dilution (Cip) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,2)));
ax6_legends{3} = strcat("1,000x Seed Dilution (Kan) k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,3)));

legend(ax6_legends,'Location','northeast','orientation','vertical');
ax6.FontSize = 14;

% Save plot
fig6_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot3normed_k2Estimation.png");
saveas(f6,fig6_plot_name);



%% Assume P_0 values and fit f_1 values

%%%%%%%%%%%%%%%%%%%%%%% Make assumptions about P_0 %%%%%%%%%%%%%%%%%%%%%%%

initial_persister_percent_1000x = [10^-1.512, 10^-2.815, 10^-5.379]; % Mean values from Fig 3d (Amp, Cip, Kan), units 1/h

ipfs_1000x = initial_persister_percent_1000x/100; % Initial persister fractions

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Guess f_1 value %%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1_initialguess =  [10^-2.708, 10^-3.691, 10^-6.052]; % Mean values from Fig 3a (Amp, Cip, Kan), units 1/h

% Initialize fitting results storage
f1_fitted = nan(1,n_cults);
f1_fit_overall_RMSE = nan(1,n_cults);
f1_fit_CI95 = nan(2,n_cults);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit f_1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation start and end times
    t_0 = 0; % hours
    t_end = 7; % hours
    t_span = [t_0, t_end];


for i = 1:n_cults

    % Save necessary info to a file that the optimization function can load
    data_passer = {t_span,t_mean_stddev{i},k1_ests(i),k2_ests(i),ipfs_1000x(i)};
    save('data_passer.mat','data_passer');

    % Run optimization
    x_guess = f1_initialguess(i);
    [x_fitted,sum_squared_residuals,residuals,exit_flag,~,~,lsq_Jacobian] = lsqnonlin(@fn_for_optimization_f1only,x_guess);

    % Calculate RMSE of fit
    RMSE_of_overall_fit = sqrt(sum_squared_residuals/length(residuals));
    f1_fit_overall_RMSE(i) = RMSE_of_overall_fit;

    % Calculate 95% CIs for fitted parameters
    CI_95 = nlparci(x_fitted,residuals,"jacobian",lsq_Jacobian);
    f1_fit_CI95(:,i) = CI_95;

    % Extract fitted parameters
    f1_fitted(i) = x_fitted;

end


%% Simulate final f_1 fitting results

% Create storage arrays for ODE solver solutions
    t_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
    log10_SPN_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
                                       % each array entry is 3 columns: S, P, and N.

% Loop through ODE solver calls
for i = 1:n_cults
    % Calculate initial conditions
    N_0 = 10^t_mean_stddev{i}(1,2); % Total cells
    P_0 = N_0*ipfs_1000x(i); % Persister cells
    S_0 = N_0-P_0; % Susceptible cells

    % Simulation
    init_conds = [S_0, P_0];
    [t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k1_ests(i),k2_ests(i),f1_fitted(i)), t_span, init_conds);

    % Store results
    t_ODEsolver{i} = t_ode;
    log10_SPN_ODEsolver{i} = log10([y_ode(:,1), y_ode(:,2), y_ode(:,1)+y_ode(:,2)]);
end


%% Plot final f_1 fitting results (Raw)

% Figure size
fig7_pixels_width=600;
fig7_pixels_height=500;

% Create figure and axis
f7=figure(7);
f7.Position = [0, 0, fig7_pixels_width, fig7_pixels_height];
ax7 = gca;

hold on;

% Plot the ODE simulation results (N)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,3),'Color',colors_RGB(i,:),'LineStyle','-','LineWidth',1);
end

% Plot the ODE simulation results (P)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,2),'Color',colors_RGB(i,:),'LineStyle','-','LineWidth',1);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',colors_RGB(i,:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(strcat("f_1 Fitting Results (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Mean log_{10}(CFU/mL)");
xlim([0,7]);
xticks(0:7);
ylim([-0.3,9]);
yticks(0:9);
grid on;

% Generate legend entries
ax7_legends = cell(n_cults,1);
ax7_legends{1} = "1,000x Seed Dilution (Amp)";
ax7_legends{2} = "1,000x Seed Dilution (Cip)";
ax7_legends{3} = "1,000x Seed Dilution (Kan)";

legend(ax7_legends,'Location','northeast','orientation','vertical');

ax7.FontSize = 14;

% Save plot
fig7_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4_FinalFit.png");
saveas(f7,fig7_plot_name);


%% Plot final f_1 fitting results (Normalized to 1)

% Figure size
fig8_pixels_width=600;
fig8_pixels_height=500;

% Create figure and axis
f8=figure(8);
f8.Position = [0, 0, fig8_pixels_width, fig8_pixels_height];
ax8 = gca;

hold on;

% Plot the ODE simulation results (N)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,3)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle','-','LineWidth',3);
end

% Plot the ODE simulation results (P)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,2)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',3);
end

% Scatter plot of mean log10 CFU/mL values
point_size = 100;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,point_size,'filled','MarkerFaceColor',colors_RGB(i,:),'MarkerEdgeColor','k','LineWidth',1,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,t_mean_stddev{i}(:,3),'color','k','LineWidth',1,'LineStyle','none','CapSize',12);
end

% title(strcat("f_1 Fitting Results (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,7]);
xticks(0:7);
ylim([-7,2.2]);
yticks(-7:1:2);
yticklabels(["10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
grid on;
set(ax8,'GridAlpha',0.4);

% Generate legend entries
ax8_legends = cell(n_cults,1);
ax8_legends{1} = "1,000x Seed Dilution (Amp)";
ax8_legends{2} = "1,000x Seed Dilution (Cip)";
ax8_legends{3} = "1,000x Seed Dilution (Kan)";

legend(ax8_legends,'Location','northeast','orientation','vertical');

ax8.FontSize = 28;

% Save plot
fig8_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4normed_FinalFit_Formatted.png");
saveas(f8,fig8_plot_name);




%% Functions

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end


function residuals = fn_for_optimization_f1only(x)

% Extract parameters from x
    f_1 = x;

% Load data from main code
    data_passer = load('data_passer.mat');

    t_span = data_passer.data_passer{1};
    t_mean_stddev = data_passer.data_passer{2};
    k_1 = data_passer.data_passer{3};
    k_2 = data_passer.data_passer{4};
    init_pers_frac_1000x = data_passer.data_passer{5};

% Run simulations and residual calculations
    residuals = [];

    % Calculate initial conditions
    N_0 = 10^t_mean_stddev(1,2); % Total cells
    P_0 = N_0*init_pers_frac_1000x; % Persister cells
    S_0 = N_0-P_0; % Susceptible cells

    % Simulation
    init_conds = [S_0, P_0];
    ode_soln = ode45(@(t,y) model_ODEs(t,y,k_1,k_2,f_1),t_span,init_conds);

    % Residual calculations
    for j = 1:length(t_mean_stddev(:,1)) 
        experimental_result = t_mean_stddev(j,2);
        simulation_state_values = deval(ode_soln,t_mean_stddev(j,1));
        simulation_result = log10(sum(simulation_state_values(1:2))); %S+P=N
        residuals = [residuals; simulation_result-experimental_result];
    end

end