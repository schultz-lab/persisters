%% Notes: Comment in either lines 4, 8, and 1293 or lines 7 and 15 to run
% the fits for all datasets or just one dataset, respectively.

for experiment_index = 1:19

%% Clear
%clear;
clearvars -except experiment_index
clc;
close all; % figures

%% Select Dataset to Analyze

% CHOOSE ONE DATAFILE TO ANALYZE
%experiment_index = 4; % Choose an integer
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
datafile = datafile_list(experiment_index);


experiment_nickname_list = ["Amp_WT_Rep1";        % 1  = Amp_WT_Rep1
                            "Amp_WT_Rep2";        % 2  = Amp_WT_Rep2
                            "Amp_WT_Rep3";        % 3  = Amp_WT_Rep3
                            "Cip_WT_Rep1";        % 4  = Cip_WT_Rep1
                            "Cip_WT_Rep2";        % 5  = Cip_WT_Rep2
                            "Cip_WT_Rep3";        % 6  = Cip_WT_Rep3
                            "Cip_WT_Rep4";        % 7  = Cip_WT_Rep4
                            "Kan_WT_Rep1";        % 8  = Kan_WT_Rep1
                            "Kan_WT_Rep2";        % 9  = Kan_WT_Rep2
                            "Kan_WT_Rep3";        % 10 = Kan_WT_Rep3
                            "Kan_WT_Rep4";        % 11 = Kan_WT_Rep4
                            "Amp_WT_Turb";        % 12 = Amp_WT_Turb
                            "Amp_ppk_Rep1";       % 13 = Amp_ppk_Rep1
                            "Amp_ppk_Rep2";       % 14 = Amp_ppk_Rep2
                            "Amp_ppk_Rep3";       % 15 = Amp_ppk_Rep3
                            "Amp_relAspoT_Rep1";  % 16 = Amp_relAspoT_Rep1
                            "Amp_relAspoT_Rep2";  % 17 = Amp_relAspoT_Rep2
                            "Amp_relAspoT_Rep3";  % 18 = Amp_relAspoT_Rep3
                            "Amp_relAspoT_Rep4"]; % 19 = Amp_relAspoT_Rep4
experiment_nickname = experiment_nickname_list(experiment_index);


%% Load information about the selected dataset

% Which antibiotic was used
drug_list = ["Amp";  % 1  = Amp_WT_Rep1
             "Amp";  % 2  = Amp_WT_Rep2
             "Amp";  % 3  = Amp_WT_Rep3
             "Cip";  % 4  = Cip_WT_Rep1
             "Cip";  % 5  = Cip_WT_Rep2
             "Cip";  % 6  = Cip_WT_Rep3
             "Cip";  % 7  = Cip_WT_Rep4
             "Kan";  % 8  = Kan_WT_Rep1
             "Kan";  % 9  = Kan_WT_Rep2
             "Kan";  % 10 = Kan_WT_Rep3
             "Kan";  % 11 = Kan_WT_Rep4
             "Amp";  % 12 = Amp_WT_Turb
             "Amp";  % 13 = Amp_ppk_Rep1
             "Amp";  % 14 = Amp_ppk_Rep2
             "Amp";  % 15 = Amp_ppk_Rep3
             "Amp";  % 16 = Amp_relAspoT_Rep1
             "Amp";  % 17 = Amp_relAspoT_Rep2
             "Amp";  % 18 = Amp_relAspoT_Rep3
             "Amp"]; % 19 = Amp_relAspoT_Rep4
drug = drug_list(experiment_index);


% What the end time (h) of the k_1 linear region is.
t_k1_end_list = [1;    % 1  = Amp_WT_Rep1
                 1;    % 2  = Amp_WT_Rep2
                 0.5;  % 3  = Amp_WT_Rep3
                 0.5;  % 4  = Cip_WT_Rep1
                 0.5;  % 5  = Cip_WT_Rep2
                 0.5;  % 6  = Cip_WT_Rep3
                 0.5;  % 7  = Cip_WT_Rep4
                 0.5;  % 8  = Kan_WT_Rep1
                 0.5;  % 9  = Kan_WT_Rep2
                 0.5;  % 10 = Kan_WT_Rep3
                 0.5;  % 11 = Kan_WT_Rep4
                 1;    % 12 = Amp_WT_Turb
                 0.5;  % 13 = Amp_ppk_Rep1
                 0.5;  % 14 = Amp_ppk_Rep2
                 0.5;  % 15 = Amp_ppk_Rep3
                 0.5;  % 16 = Amp_relAspoT_Rep1
                 1;    % 17 = Amp_relAspoT_Rep2
                 0.5;  % 18 = Amp_relAspoT_Rep3
                 0.5]; % 19 = Amp_relAspoT_Rep4
t_k1_end = t_k1_end_list(experiment_index);
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% What the start time (h) of the k_2 linear region is
t_k2_start_list = [3;    % 1  = Amp_WT_Rep1
                   2;    % 2  = Amp_WT_Rep2
                   1.5;  % 3  = Amp_WT_Rep3
                   3;    % 4  = Cip_WT_Rep1
                   2;    % 5  = Cip_WT_Rep2
                   2;    % 6  = Cip_WT_Rep3
                   2;    % 7  = Cip_WT_Rep4
                   1.5;  % 8  = Kan_WT_Rep1
                   2;    % 9  = Kan_WT_Rep2
                   2;    % 10 = Kan_WT_Rep3
                   2;    % 11 = Kan_WT_Rep4
                   2;    % 12 = Amp_WT_Turb
                   1.5;  % 13 = Amp_ppk_Rep1
                   1.5;  % 14 = Amp_ppk_Rep2
                   2;    % 15 = Amp_ppk_Rep3
                   2;    % 16 = Amp_relAspoT_Rep1
                   2;    % 17 = Amp_relAspoT_Rep2
                   2;    % 18 = Amp_relAspoT_Rep3
                   0.5]; % 19 = Amp_relAspoT_Rep4
t_k2_start = t_k2_start_list(experiment_index);
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% Set of dilution factors used
dil_fact_set_list = [100, 1000, 10000;    % 1  = Amp_WT_Rep1
                     100, 1000, 10000;    % 2  = Amp_WT_Rep2
                     100, 1000, 10000;    % 3  = Amp_WT_Rep3
                     100, 1000, 10000;    % 4  = Cip_WT_Rep1
                     100, 1000, 10000;    % 5  = Cip_WT_Rep2
                     100, 1000, 10000;    % 6  = Cip_WT_Rep3
                     100, 1000, 10000;    % 7  = Cip_WT_Rep4
                     100, 1000, 10000;    % 8  = Kan_WT_Rep1
                     100, 1000, 10000;    % 9  = Kan_WT_Rep2
                     100, 1000, 10000;    % 10 = Kan_WT_Rep3
                     100, 1000, 10000;    % 11 = Kan_WT_Rep4
                      20,  200,  2000;    % 12 = Amp_WT_Turb
                     100, 1000, 10000;    % 13 = Amp_ppk_Rep1
                     100, 1000, 10000;    % 14 = Amp_ppk_Rep2
                     100, 1000, 10000;    % 15 = Amp_ppk_Rep3
                     100, 1000, 10000;    % 16 = Amp_relAspoT_Rep1
                     100, 1000, 10000;    % 17 = Amp_relAspoT_Rep2
                     100, 1000, 10000;    % 18 = Amp_relAspoT_Rep3
                     100, 1000, 10000];   % 19 = Amp_relAspoT_Rep4
dil_fact_set = dil_fact_set_list(experiment_index,:);


% First Culture Type
first_cult_type_list = ["batch";    % 1  = Amp_WT_Rep1
                        "batch";    % 2  = Amp_WT_Rep2
                        "batch";    % 3  = Amp_WT_Rep3
                        "batch";    % 4  = Cip_WT_Rep1
                        "batch";    % 5  = Cip_WT_Rep2
                        "batch";    % 6  = Cip_WT_Rep3
                        "batch";    % 7  = Cip_WT_Rep4
                        "batch";    % 8  = Kan_WT_Rep1
                        "batch";    % 9  = Kan_WT_Rep2
                        "batch";    % 10 = Kan_WT_Rep3
                        "batch";    % 11 = Kan_WT_Rep4
                        "turbidostat";   % 12 = Amp_WT_Turb
                        "batch";    % 13 = Amp_ppk_Rep1
                        "batch";    % 14 = Amp_ppk_Rep2
                        "batch";    % 15 = Amp_ppk_Rep3
                        "batch";    % 16 = Amp_relAspoT_Rep1
                        "batch";    % 17 = Amp_relAspoT_Rep2
                        "batch";    % 18 = Amp_relAspoT_Rep3
                        "batch"];   % 19 = Amp_relAspoT_Rep4
first_cult_type = first_cult_type_list(experiment_index);


% The genotype
genotype_list = ["WT";    % 1  = Amp_WT_Rep1
                 "WT";    % 2  = Amp_WT_Rep2
                 "WT";    % 3  = Amp_WT_Rep3
                 "WT";    % 4  = Cip_WT_Rep1
                 "WT";    % 5  = Cip_WT_Rep2
                 "WT";    % 6  = Cip_WT_Rep3
                 "WT";    % 7  = Cip_WT_Rep4
                 "WT";    % 8  = Kan_WT_Rep1
                 "WT";    % 9  = Kan_WT_Rep2
                 "WT";    % 10 = Kan_WT_Rep3
                 "WT";    % 11 = Kan_WT_Rep4
                 "WT";    % 12 = Amp_WT_Turb
                 "ppk";    % 13 = Amp_ppk_Rep1
                 "ppk";    % 14 = Amp_ppk_Rep2
                 "ppk";    % 15 = Amp_ppk_Rep3
                 "relAspoT";    % 16 = Amp_relAspoT_Rep1
                 "relAspoT";    % 17 = Amp_relAspoT_Rep2
                 "relAspoT";    % 18 = Amp_relAspoT_Rep3
                 "relAspoT"];   % 19 = Amp_relAspoT_Rep4
genotype = genotype_list(experiment_index);



%% Define a color palette for plots
colors_RGB = [1.0, 0.0, 0.0; % red
              0.3, 0.7, 0.3; % green
              0.0, 0.0, 1.0]; % blue


%% Load Dataset & Extract Data

% Load data
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
% each dilution condition
t_mean_stddev = {};
for i = 1:n_cults
    t_mean_stddev{i} = data(1:last_time_index(i),[1,1+i,1+n_cults+i]);
end

% Skip over missing data in particular datasets
if experiment_nickname=="Amp_relAspoT_Rep1"
    t_mean_stddev{2} = data([1:3,5:9],[1,1+2,1+n_cults+2]);
    t_mean_stddev{3} = data([1:3,5:8],[1,1+3,1+n_cults+3]);
end


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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax1_legends = cell(n_cults,1);
for i = 1:n_cults
    ax1_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

% Generate legend entries
ax2_legends = cell(n_cults,1);
for i = 1:n_cults
    ax2_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

legend(ax2_legends,'Location','northeast','orientation','vertical');
ax2.FontSize = 14;

% Save plot
fig2_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot1normed_ExperimentalData.png");
saveas(f2,fig2_plot_name);



%% Estimate k_1 from the first linear region

% Truncate data to the first linear region
t_mean_stddev_k1est = cell(size(t_mean_stddev));
for i = 1:n_cults
    n_tpts = sum(t_mean_stddev{i}(:,1) <= t_k1_end);
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

% Calculate the mean estimated k_1 value and its stddev
k1_est_mean = mean(k1_fits(3,:));
k1_est_stddev = std(k1_fits(3,:));


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
    plot(0:0.01:9,polyval(linear_k1fit_coefs_storage(i,:),0:0.01:9),'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax3_legends = cell(n_cults,1);
for i = 1:n_cults
    ax3_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_fits(2,i)));
end

legend(ax3_legends,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{1}: %.3f ± %.3f h^{-1}",k1_est_mean,k1_est_stddev),"FontSize",14);
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
    plot(0:0.01:9,polyval(linear_k1fit_coefs_storage(i,:),0:0.01:9)-t_mean_stddev{i}(1,2)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

% Generate legend entries
ax4_legends = cell(n_cults,1);
for i = 1:n_cults
    ax4_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution k_1 Fit Line: RMSE=",sprintf("%.3f log_{10}(%% Survival)",k1_fits(2,i)));
end

legend(ax4_legends,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{1}: %.3f ± %.3f h^{-1}",k1_est_mean,k1_est_stddev),"FontSize",14);
ax4.FontSize = 14;

% Save plot
fig4_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot2normed_k1Estimation.png");
saveas(f4,fig4_plot_name);



%% Estimate k_2 from the second linear region

% Truncate data to the second linear region
t_mean_stddev_k2est = cell(size(t_mean_stddev));
for i = 1:n_cults
    n_tpts = sum(t_mean_stddev{i}(:,1) >= t_k2_start);
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
k2_est_mean = mean(k2_fits(3,:));
k2_est_stddev = std(k2_fits(3,:));


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
    plot(0:0.01:9,polyval(linear_k2fit_coefs_storage(i,:),0:0.01:9),'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax5_legends = cell(n_cults,1);
for i = 1:n_cults
    ax5_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_fits(2,i)));
end

legend(ax5_legends,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{2}: %.3f ± %.3f h^{-1}",k2_est_mean,k2_est_stddev),"FontSize",14);
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
    plot(0:0.01:9,polyval(linear_k2fit_coefs_storage(i,:),0:0.01:9)-t_mean_stddev{i}(1,2)+2,'Color',colors_RGB(i,:),'LineStyle','--','LineWidth',1);
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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

% Generate legend entries
ax6_legends = cell(n_cults,1);
for i = 1:n_cults
    ax6_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution k_2 Fit Line: RMSE=",sprintf("%.3f log_{10}(%% Survival)",k2_fits(2,i)));
end

legend(ax6_legends,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{2}: %.3f ± %.3f h^{-1}",k2_est_mean,k2_est_stddev),"FontSize",14);
ax6.FontSize = 14;

% Save plot
fig6_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot3normed_k2Estimation.png");
saveas(f6,fig6_plot_name);



%% Fitting f_1 (and P_0 for batch)


if first_cult_type=="batch"
%% Simultaneously fit f_1 and P_0 values

%%%%%%%%%%%% Estimate initial_persister_fraction_leastdilute %%%%%%%%%%%%%%

% Assume that a negligible fraction of the persisters in the least dilute
% culture form as a result of the antibiotic (i.e., almost all persisters
% in the least dilute culture existed at time 0). Therefore, the intercept
% between the k_2 fit line and the vertical axis provides an estimate of
% P_0 for the least dilute culture.

% The least dilute culture has index i=1 in t_mean_stddev{i} and
% dil_fact_set(i)

% Calculate the P_0 estimate based on the mean k_2 estimate, not just the
% k_2 value from the least dilute culture's linear fit. This should make
% the fitting less sensitive to outliers in the least dilute culture's cell
% count measurements.

% Loop to find the vertical intercept (to 0.01 log10 precision) that 
% gives least squares (in (log10)^2 units) fit to the second linear region.
y_int_test_lb = round(linear_k2fit_coefs_storage(1,2)-2,2); % 2 log10 units below indiviudal linear fit
y_int_test_ub = round(linear_k2fit_coefs_storage(1,2)+2,2); % 2 log10 units below indiviudal linear fit
y_ints_testset = y_int_test_lb:0.01:y_int_test_ub;
n_y_ints = length(y_ints_testset);
SSE = zeros(size(y_ints_testset));
for i=1:n_y_ints
    % y   =         m               *           t                +        b
    y_fit = (-1*k2_est_mean/log(10))*t_mean_stddev_k2est{1}(:,1) + y_ints_testset(i);
    SSE(i) = sum((y_fit-t_mean_stddev_k2est{1}(:,2)).^2);
end
[~, y_int_index] = min(SSE);
if (y_int_index==1 || y_int_index==n_y_ints)
    disp("Warning: Initial P_0 guess poorly fitted.");
end
P_0_leastdilute_est = 10^y_ints_testset(y_int_index);
init_pers_frac_leastdilute_est = P_0_leastdilute_est/(10.^t_mean_stddev{1}(1,2));

% Then, assume that persister formation was negligigble in the "pre-drug  
% cultures" (the cultures that are grown up to a consistent OD after the 
% dilution of the overnight stationary phase culture but before the 
% addition of antibiotic). Therefore, we can assume that the persister
% fractions in the more dilute cultures are related to the persister
% fraction in the least dilute culture by the dilution factors
    % Form: init_pers_frac_moredilute_est = 
    % (dil_fact_leastdilute/dil_fact_moredilute)*init_pers_frac_leastdilute_est


%%%%%%%%%%%%%%%%%%%%%%%%% Initial guess f_1 value %%%%%%%%%%%%%%%%%%%%%%%%%
if drug=="Kan"
    f1_initialguess = 10^-7; % units: 1/h
else
    f1_initialguess = 10^-3; % units: 1/h
end

%%%%%%%%% Simultaneously fit f_1 and initial persister fraction %%%%%%%%%%%

% Set simulation start and end times
    t_0 = 0; % hours
    if experiment_nickname=="Amp_WT_Rep3"
        t_end = 9; % hours
    else
        t_end = 8; % hours
    end
    t_span = [t_0, t_end];

% Save necessary info to a file that the optimization function can load
    data_passer = {t_span,t_mean_stddev,k1_est_mean,k2_est_mean,dil_fact_set};
    save('data_passer.mat','data_passer');

% Set optimization bounds
    f1_lb = -inf;
    f1_ub = inf;
    ipf_ld_lb = -inf;
    ipf_ld_ub = inf;

    x_lb = [f1_lb; ipf_ld_lb];
    x_ub = [f1_ub; ipf_ld_ub];

% Run optimization
    x_guess = [f1_initialguess; init_pers_frac_leastdilute_est];
    [x_fitted,sum_squared_residuals,residuals,exit_flag,~,~,lsq_Jacobian] = lsqnonlin(@fn_for_optimization_P0andf1,x_guess,x_lb,x_ub);

% Calculate RMSE of fit
    RMSE_of_overall_fit = sqrt(sum_squared_residuals/length(residuals));

% Calculate 95% CIs for fitted parameters
    CI_95 = nlparci(x_fitted,residuals,"jacobian",lsq_Jacobian);

% Extract fitted parameters
    f1_fitted = x_fitted(1);
    init_pers_frac_leastdilute_fitted = x_fitted(2);



%% Simulate final f_1 and P_0 fitting results 

% Create storage arrays for ODE solver solutions
    t_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
    log10_SPN_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
                                       % each array entry is 3 columns: S, P, and N.

% Loop through ODE solver calls
for i = 1:n_cults
    % Calculate initial conditions
    N_0 = 10^t_mean_stddev{i}(1,2); % Total cells
    P_0 = N_0*(dil_fact_set(1)/dil_fact_set(i))*init_pers_frac_leastdilute_fitted; % Total persister cells
    S_0 = N_0-P_0; % Susceptible cells

    % Simulation
    init_conds = [S_0, P_0];
    [t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k1_est_mean,k2_est_mean,f1_fitted), t_span, init_conds);
    
    % Store results
    t_ODEsolver{i} = t_ode;
    log10_SPN_ODEsolver{i} = log10([y_ode(:,1), y_ode(:,2), y_ode(:,1)+y_ode(:,2)]);
end


%% Plot final f_1 and P_0 fitting results (Raw)

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

title(strcat("f_1 and P_0 Simulataneous Fitting Results (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Mean log_{10}(CFU/mL)");
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax7_legends = cell(n_cults,1);
for i = 1:n_cults
    ax7_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

legend(ax7_legends,'Location','northeast','orientation','vertical');
text(2.5,7.25,sprintf("Fitted: f_{1} = %.2e h^{-1}",f1_fitted),"FontSize",14,"FontWeight","bold");
text(2.5,6.75,sprintf("95%% CI: [%.2e, %.2e] h^{-1}",CI_95(1,:)),"FontSize",14);
if init_pers_frac_leastdilute_fitted>0
    text(2.5,6,sprintf("Fitted: P_{0,%dx} = %.2f log_{10}(CFU/mL)",dil_fact_set(1),log10((10^t_mean_stddev{1}(1,2))*init_pers_frac_leastdilute_fitted)),"FontSize",14,"FontWeight","bold");
else
    text(2.5,6,sprintf("Fitted: P_{0,%dx} = NEG DIR %.2f log_{10}(CFU/mL)",dil_fact_set(1),log10((10^t_mean_stddev{1}(1,2))*abs(init_pers_frac_leastdilute_fitted))),"FontSize",14,"FontWeight","bold");
end
if CI_95(2,1)>0 && CI_95(2,2)>0
    text(2.5,5.5,sprintf("95%% CI: [%.2f, %.2f] log_{10}(CFU/mL)",log10((10^t_mean_stddev{1}(1,2))*(CI_95(2,:)))),"FontSize",14);
elseif CI_95(2,1)<=0 && CI_95(2,2)>0
    text(2.5,5.5,sprintf("95%% CI: [NEG DIR %.2f, %.2f] log_{10}(CFU/mL)",log10((10^t_mean_stddev{1}(1,2))*abs((CI_95(2,:))))),"FontSize",14);
elseif CI_95(2,1)<=0 && CI_95(2,2)<=0
    text(2.5,5.5,sprintf("95%% CI: [NEG DIR %.2f, NEG DIR %.2f] log_{10}(CFU/mL)",log10((10^t_mean_stddev{1}(1,2))*abs((CI_95(2,:))))),"FontSize",14);
end
text(2.5,4.75,sprintf("Overall RMSE = %.3f log_{10}(CFU/mL)",RMSE_of_overall_fit),"FontSize",14);

ax7.FontSize = 14;

% Save plot
fig7_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4_FinalFit.png");
saveas(f7,fig7_plot_name);


%% Plot final f_1 and P_0 fitting results (Normalized to 1)

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
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,point_size,'filled','MarkerEdgeColor',colors_RGB(i,:),'MarkerFaceColor',colors_RGB(i,:),'LineWidth',2,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{i}(1,2)+2,t_mean_stddev{i}(:,3),'color','k','LineWidth',1,'LineStyle','none','CapSize',12);
end

%title(strcat("f_1 and P_0 Simulataneous Fitting Results (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
    text_shift_vert = -1;
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
    text_shift_vert = 0;
end
grid on;
set(ax8,'GridAlpha',0.4);

% Generate legend entries
ax8_legends = cell(n_cults,1);
for i = 1:n_cults
    ax8_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

legend(ax8_legends,'Location','northeast','orientation','vertical');
ax8.FontSize = 28;

% Save plot
fig8_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4normed_FinalFit_Formatted.png");
saveas(f8,fig8_plot_name);




elseif first_cult_type=="turbidostat"
%% Assume P_0 value and fit f_1 value

%%%%%%%%%%%%%%%%%%%%%%% Make assumptions about P_0 %%%%%%%%%%%%%%%%%%%%%%%
% Assume that any persisters formed during the exponential phase
% turbidostat culture are diluted away due to their slow/zero growth rate.
% Then, if no persisters form during the "pre-drug cultures" (the cultures 
% that are grown up to a consistent OD after the dilution of the turbidostat 
% exponential phase culture into batch cultures but before the addition of 
% antibiotic), the number of persisters present when the drug is applied
% at time zero is approximately 0. This suggests a minimum assumed 
% P_0 value of 0. The minimum assumed P_0 value leads to the maximum
% fitted f_1 value.

% However, we acknowledge that a small fraction of the cells in the pre- 
% drug batch culture likely become "spontaneous persisters" (formed
% without the presence of stressors such as nutrient depletion or
% antibiotic treatment). Based on our simulation of Balaban's 2004 model in
% the Supplementary Figure 1, we assume that the initial persister
% fraction for the turbidostat-origin cultures is approximately 8*10^-7
% CFU/mL.

% By assuming the higher value of P_0, we guarantee that the fitted value
% of f_1 is lower. We find that such a "lowered" fitted f_1 value is
% still substantial in magnitude and positive in sign, providing evidence
% that drug-induced persistence is a real phenomenon.

% Note: we also assume that P_0 is the same for all dilutions of the
% turbidostat culture because we assume that there are no persisters in the
% turbidostat and the spontaneous persister fraction in the pre-drug
% culture equilibrates before the OD threshold for drug addition is met 
% (see converging purple lines in Supplementary Fig 1).

if drug == "Amp"
    init_persister_frac = 8*10^-7;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Initial guess f_1 value %%%%%%%%%%%%%%%%%%%%%%%%%
if drug=="Kan"
    f1_initialguess = 10^-7; % units: 1/h
else
    f1_initialguess = 10^-3; % units: 1/h
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fit f_1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation start and end times
    t_0 = 0; % hours
    if experiment_nickname=="Amp_WT_Rep3"
        t_end = 9; % hours
    else
        t_end = 8; % hours
    end
    t_span = [t_0, t_end];

% Save necessary info to a file that the optimization function can load
    data_passer = {t_span,t_mean_stddev,k1_est_mean,k2_est_mean,dil_fact_set,init_persister_frac};
    save('data_passer.mat','data_passer');

% Set optimization bounds
    f1_lb = -inf;
    f1_ub = inf;
    
    x_lb = f1_lb;
    x_ub = f1_ub;

% Run optimization
    x_guess = f1_initialguess;
    [x_fitted,sum_squared_residuals,residuals,exit_flag,~,~,lsq_Jacobian] = lsqnonlin(@fn_for_optimization_turbidostatf1only,x_guess,x_lb,x_ub);

% Calculate RMSE of fit
    RMSE_of_overall_fit = sqrt(sum_squared_residuals/length(residuals));

% Calculate 95% CIs for fitted parameters
    CI_95 = nlparci(x_fitted,residuals,"jacobian",lsq_Jacobian);

% Extract fitted parameters
    f1_fitted = x_fitted;



%% Simulate final f_1 fitting results

% Create storage arrays for ODE solver solutions
    t_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
    log10_SPN_ODEsolver = cell(n_cults,1); % indices increasing from least dilute to most dilute culture
                                       % each array entry is 3 columns: S, P, and N.

% Loop through ODE solver calls
for i = 1:n_cults
    % Calculate initial conditions
    N_0 = 10^t_mean_stddev{i}(1,2); % Total cells
    P_0 = N_0*init_persister_frac; % Total persister cells
    S_0 = N_0-P_0; % Susceptible cells

    % Simulation
    init_conds = [S_0, P_0];
    [t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k1_est_mean,k2_est_mean,f1_fitted), t_span, init_conds);
    
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
if experiment_nickname=="Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax7_legends = cell(n_cults,1);
for i = 1:n_cults
    ax7_legends{i} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

legend(ax7_legends,'Location','northeast','orientation','vertical');
text(2.5,7.25,sprintf("Fitted: f_{1} = %.2e h^{-1}",f1_fitted),"FontSize",14,"FontWeight","bold");
text(2.5,6.75,sprintf("95%% CI: [%.2e, %.2e] h^{-1}",CI_95(1,:)),"FontSize",14);
text(2.5,6,sprintf("Overall RMSE = %.3f log_{10}(CFU/mL)",RMSE_of_overall_fit),"FontSize",14);

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

linewidths = [3.6,3,2.4];
% Plot the ODE simulation results (P)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,2)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle',':','LineWidth',linewidths(i));
end

linewidths = [7,5,1.5];
% Plot the ODE simulation results (N)
for i = 1:n_cults
    plot(t_ODEsolver{i},log10_SPN_ODEsolver{i}(:,3)-log10_SPN_ODEsolver{i}(1,3)+2,'Color',colors_RGB(i,:),'LineStyle','-','LineWidth',linewidths(i));
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

%title(strcat("f_1 Fitting Results (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
if experiment_nickname == "Amp_WT_Rep3"
    xlim([0,9]);
    xticks(0:9);
else
    xlim([0,8]);
    xticks(0:8);
end
if drug=="Kan"
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
    text_shift_vert = -1;
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
    text_shift_vert = 0;
end
grid on;
set(ax8,'GridAlpha',0.4);

% Generate legend entries
ax8_legends = cell(n_cults,1);
for i = 1:n_cults
    ax8_legends{i} = "";
    ax8_legends{i+n_cults} = strcat(num2str(dil_fact_set(i)),"x Seed Dilution");
end

legend(ax8_legends,'Location','northeast','orientation','vertical','FontSize',18);

ax8.FontSize = 28;

% Save plot
fig8_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4normed_FinalFit_Formatted.png");
saveas(f8,fig8_plot_name);



end




%% Save key results to a table

% Columns:
% experiment_index
% experiment_nickname
% drug
% first_cult_type
% dil_fact_set
% genotype
% t_k1_end
% k1_est_mean
% k1_est_stddev
% t_k2_start
% k2_est_mean
% k2_est_stddev
% f1_fitted
% f1_95CI_lb = CI_95(1,1);
% f1_95CI_ub = CI_95(1,2);
% init_pers_frac_leastdilute_fitted
% ipf_leastdilute_95CI_lb = CI_95(2,1);
% ipf_leastdilute_95CI_ub = CI_95(2,2);
% RMSE_of_overall_fit


f1_95CI_lb = CI_95(1,1);
f1_95CI_ub = CI_95(1,2);


if (first_cult_type=="turbidostat")
    init_pers_frac_leastdilute_fitted = NaN;
    ipf_leastdilute_95CI_lb = NaN;
    ipf_leastdilute_95CI_ub = NaN;
else
    ipf_leastdilute_95CI_lb = CI_95(2,1);
    ipf_leastdilute_95CI_ub = CI_95(2,2);
end


fit_results_tab = table(experiment_index, experiment_nickname, drug,...
    first_cult_type, dil_fact_set, genotype, t_k1_end, k1_est_mean,...
    k1_est_stddev, t_k2_start, k2_est_mean, k2_est_stddev, f1_fitted,...
    f1_95CI_lb, f1_95CI_ub, init_pers_frac_leastdilute_fitted,...
    ipf_leastdilute_95CI_lb, ipf_leastdilute_95CI_ub, RMSE_of_overall_fit);

save(strcat("Output_Tables/Experiment",num2str(experiment_index),"_ResultsTable.mat"),"fit_results_tab");



end



%% Functions

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end


function residuals = fn_for_optimization_P0andf1(x)

% Extract parameters from x
    f_1 = x(1);
    init_pers_frac_leastdilute = x(2);

% Load data from main code
    data_passer = load('data_passer.mat');

    t_span = data_passer.data_passer{1};
    t_mean_stddev = data_passer.data_passer{2};
    k_1 = data_passer.data_passer{3};
    k_2 = data_passer.data_passer{4};
    dil_fact_set = data_passer.data_passer{5};

% Loop through simulations and residual calculations
    residuals = [];
    for i=1:length(t_mean_stddev)
        % Calculate initial conditions
        N_0 = 10^t_mean_stddev{i}(1,2); % Total cells
        P_0 = N_0*(dil_fact_set(1)/dil_fact_set(i))*init_pers_frac_leastdilute; % Persister cells
        S_0 = N_0-P_0; % Susceptible cells

        % Simulation
        init_conds = [S_0, P_0];
        ode_soln = ode45(@(t,y) model_ODEs(t,y,k_1,k_2,f_1),t_span,init_conds);
        
        % Residual calculations
        for j = 1:length(t_mean_stddev{i}(:,1))
            experimental_result = t_mean_stddev{i}(j,2);
            simulation_state_values = deval(ode_soln,t_mean_stddev{i}(j,1));
            simulation_result = log10(sum(simulation_state_values(1:2))); %S+P=N
            residuals = [residuals; simulation_result-experimental_result];
        end
    end
end


function residuals = fn_for_optimization_turbidostatf1only(x)

% Extract parameters from x
    f_1 = x(1);

% Load data from main code
    data_passer = load('data_passer.mat');

    t_span = data_passer.data_passer{1};
    t_mean_stddev = data_passer.data_passer{2};
    k_1 = data_passer.data_passer{3};
    k_2 = data_passer.data_passer{4};
    dil_fact_set = data_passer.data_passer{5};
    init_pers_frac = data_passer.data_passer{6};

% Loop through simulations and residual calculations
    residuals = [];
    for i=1:length(t_mean_stddev)
        % Calculate initial conditions
        N_0 = 10^t_mean_stddev{i}(1,2); % Total cells
        P_0 = N_0*init_pers_frac; % Persister cells
        S_0 = N_0-P_0; % Susceptible cells

        % Simulation
        init_conds = [S_0, P_0];
        ode_soln = ode45(@(t,y) model_ODEs(t,y,k_1,k_2,f_1),t_span,init_conds);
        
        % Residual calculations
        for j = 1:length(t_mean_stddev{i}(:,1))
            experimental_result = t_mean_stddev{i}(j,2);
            simulation_state_values = deval(ode_soln,t_mean_stddev{i}(j,1));
            simulation_result = log10(sum(simulation_state_values(1:2))); %S+P=N
            residuals = [residuals; simulation_result-experimental_result];
        end
    end
end