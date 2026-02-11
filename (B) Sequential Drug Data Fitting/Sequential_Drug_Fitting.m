%% Notes
% This code assumes the data for the first drug applied is in the first
% columns of mean and stddev data in the "{whatever}_Data.txt" file.

% Comment in either lines 18, 22, and 1151 or lines 21 and 29 to run
% the fits for all datasets or just one dataset, respectively.

% This code was written to accomodate only 1 second drug, so datasets with
% 2 second drugs (an antibiotic and a PBS control) have been split into 
% pairs of two datasets. The drug in parentheses in the experiment nickname
% is the drug NOT included in that split dataset.

% Lines 1078-1113 have save() statements that can be commented in and 
% renamed as needed to save the final plotted results. These saved results
% are used to make the "Combined" images with both second drugs (antibiotic
% and PBS) in the  "Combined_..." folders.

for experiment_index = 1:12

%% Clear
% clear;
clearvars -except experiment_index
clc;
close all; % figures

%% Select Dataset to Analyze

% CHOOSE ONE DATAFILE TO ANALYZE
% experiment_index = 5; % Choose an integer
                                                         % Experiment
                                                         % Indices
datafile_list = ["Seqdrug_AmpKan_081725_Data.txt";       %  1 = AmpKan
                 "Seqdrug_KanAmp(PBS)_100525_Data.txt";  %  2 = KanAmp(PBS)
                 "Seqdrug_Kan(Amp)PBS_100525_Data.txt";  %  3 = Kan(Amp)PBS
                 "Seqdrug_AmpCip_122325_Data.txt";       %  4 = AmpCip
                 "Seqdrug_CipAmp(PBS)_122225_Data.txt";  %  5 = CipAmp(PBS)
                 "Seqdrug_Cip(Amp)PBS_122225_Data.txt";  %  6 = Cip(Amp)PBS
                 "Seqdrug_CipKan_122225_Data.txt";       %  7 = CipKan
                 "Seqdrug_KanCip(PBS)_121625_Data.txt";  %  8 = KanCip(PBS)
                 "Seqdrug_Kan(Cip)PBS_121625_Data.txt";  %  9 = Kan(Cip)PBS
                 "Seqdrug_AmpCip1000x_110725_Data.txt";  % 10 = AmpCip1000x
                 "Seqdrug_CipAmp1000x_120525_Data.txt";  % 11 = CipAmp1000x
                 "Seqdrug_CipKan1000x_120525_Data.txt"]; % 12 = CipKan1000x
datafile = datafile_list(experiment_index);


experiment_nickname_list = ["AmpKan";       %  1 = AmpKan
                            "KanAmp(PBS)";  %  2 = KanAmp(PBS)
                            "Kan(Amp)PBS";  %  3 = Kan(Amp)PBS
                            "AmpCip";       %  4 = AmpCip
                            "CipAmp(PBS)";  %  5 = CipAmp(PBS)
                            "Cip(Amp)PBS";  %  6 = Cip(Amp)PBS
                            "CipKan";       %  7 = CipKan
                            "KanCip(PBS)";  %  8 = KanCip(PBS)
                            "Kan(Cip)PBS";  %  9 = Kan(Cip)PBS
                            "AmpCip1000x";  % 10 = AmpCip1000x
                            "CipAmp1000x";  % 11 = CipAmp1000x
                            "CipKan1000x"]; % 12 = CipKan1000x
experiment_nickname = experiment_nickname_list(experiment_index);


%% Load information about the selected dataset

% Which antibiotics were used in what order
drugs_list = ["Amp", "Kan";  %  1 = AmpKan
              "Kan", "Amp";  %  2 = KanAmp(PBS)
              "Kan", "PBS";  %  3 = Kan(Amp)PBS
              "Amp", "Cip";  %  4 = AmpCip
              "Cip", "Amp";  %  5 = CipAmp(PBS)
              "Cip", "PBS";  %  6 = Cip(Amp)PBS
              "Cip", "Kan";  %  7 = CipKan
              "Kan", "Cip";  %  8 = KanCip(PBS)
              "Kan", "PBS";  %  9 = Kan(Cip)PBS  
              "Amp", "Cip";  % 10 = AmpCip1000x
              "Cip", "Amp";  % 11 = CipAmp1000x
              "Cip", "Kan"]; % 12 = CipKan1000x
drugs = drugs_list(experiment_index,:);


% What the end time (h) of the first drug's k_1 linear region is.
t_k1_drug1_end_list = [0.5;   %  1 = AmpKan
                       0.5;   %  2 = KanAmp(PBS)
                       0.5;   %  3 = Kan(Amp)PBS % Must be same as KanAmp(PBS) (same data)
                       0.5;   %  4 = AmpCip
                       0.5;   %  5 = CipAmp(PBS) 
                       0.5;   %  6 = Cip(Amp)PBS % Must be same as CipAmp(PBS) (same data)
                       0.5;   %  7 = CipKan % Must be same as CipAmp(PBS) (same data)
                       0.5;   %  8 = KanCip(PBS)
                       0.5;   %  9 = Kan(Cip)PBS % Must be same as KanCip(PBS) (same data)
                       0.5;   % 10 = AmpCip1000x
                       0.5;   % 11 = CipAmp1000x
                       0.5];  % 12 = CipKan1000x % Must be same as CipAmp1000x (same data)
t_k1_drug1_end = t_k1_drug1_end_list(experiment_index);
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% What the start time (h) of the first drug's k_2 linear region is
t_k2_drug1_start_list = [2;   %  1 = AmpKan
                         1;   %  2 = KanAmp(PBS)
                         1;   %  3 = Kan(Amp)PBS % Must be same as KanAmp(PBS) (same data)
                         2;   %  4 = AmpCip
                         1;   %  5 = CipAmp(PBS) 
                         1;   %  6 = Cip(Amp)PBS % Must be same as CipAmp(PBS) (same data)
                         1;   %  7 = CipKan % Must be same as CipAmp(PBS) (same data)
                         1;   %  8 = KanCip(PBS)
                         1;   %  9 = Kan(Cip)PBS % Must be same as KanCip(PBS) (same data)
                         2;   % 10 = AmpCip1000x
                         1;   % 11 = CipAmp1000x
                         1];  % 12 = CipKan1000x % Must be same as CipAmp1000x (same data)
t_k2_drug1_start = t_k2_drug1_start_list(experiment_index);
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% What the end time (h) of the second drug's k_1 linear region is.
t_k1_drug2_end_list = [3;   %  1 = AmpKan
                       inf; %  2 = KanAmp(PBS)
                       inf; %  3 = Kan(Amp)PBS
                       3;   %  4 = AmpCip
                       inf; %  5 = CipAmp(PBS)
                       inf; %  6 = Cip(Amp)PBS
                       3;   %  7 = CipKan
                       inf; %  8 = KanCip(PBS)
                       inf; %  9 = Kan(Cip)PBS
                       3;   % 10 = AmpCip1000x
                       inf; % 11 = CipAmp1000x
                       4];  % 12 = CipKan1000x
t_k1_drug2_end = t_k1_drug2_end_list(experiment_index);
% NOTE: Values of "inf" are used when there is assumed to be only a k_2
% killing phase for the second drug. That is to say, it is assumed that all
% cells remaining at the time of the second drug's addition are persistent
% to the second drug. In this case, the code doesn't actually group all
% remaining cells into a pesister group. Instead, the susceptible cells and
% persister cells are calculated normally and the susceptible cells' k_1 is
% set to the k_2 value. This is okay because we assume growth and f_1 are 
% zero in these sequential drug fits (f_1=0 assumption commented on further 
% below in the code section "Fitting Curves to Experimental Data").
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% What the start time (h) of the second drug's k_2 linear region is
t_k2_drug2_start_list = [4;   %  1 = AmpKan
                         2;   %  2 = KanAmp(PBS)
                         2;   %  3 = Kan(Amp)PBS
                         5;   %  4 = AmpCip
                         2;   %  5 = CipAmp(PBS)
                         2;   %  6 = Cip(Amp)PBS
                         4;   %  7 = CipKan
                         2;   %  8 = KanCip(PBS)
                         2;   %  9 = Kan(Cip)PBS
                         3;   % 10 = AmpCip1000x
                         2;   % 11 = CipAmp1000x
                         4];  % 12 = CipKan1000x
t_k2_drug2_start = t_k2_drug2_start_list(experiment_index);
% Note: Timepoint selected for each experiment to qualitatively balance 3 
% factors: (1) Achieving low RMSE values for each of the linear fits, 
% (2) inclusion of as many datapoints as possible in the linear region, 
% and (3) achieving a low RMSE for the overall (f1 and P0) fit


% Set of dilution factors used
dil_fact_set_list = [100;   %  1 = AmpKan
                     100;   %  2 = KanAmp(PBS)
                     100;   %  3 = Kan(Amp)PBS
                     100;   %  4 = AmpCip
                     100;   %  5 = CipAmp(PBS)
                     100;   %  6 = Cip(Amp)PBS
                     100;   %  7 = CipKan
                     100;   %  8 = KanCip(PBS)
                     100;   %  9 = Kan(Cip)PBS
                     1000;  % 10 = AmpCip1000x
                     1000;  % 11 = CipAmp1000x
                     1000]; % 12 = CipKan1000x
dil_fact_set = dil_fact_set_list(experiment_index,:);


% First Culture Type
first_cult_type_list = ["batch";   %  1 = AmpKan
                        "batch";   %  2 = KanAmp(PBS)
                        "batch";   %  3 = Kan(Amp)PBS
                        "batch";   %  4 = AmpCip
                        "batch";   %  5 = CipAmp(PBS)
                        "batch";   %  6 = Cip(Amp)PBS
                        "batch";   %  7 = CipKan
                        "batch";   %  8 = KanCip(PBS)
                        "batch";   %  9 = Kan(Cip)PBS
                        "batch";   % 10 = AmpCip1000x
                        "batch";   % 11 = CipAmp1000x
                        "batch"];  % 12 = CipKan1000x
first_cult_type = first_cult_type_list(experiment_index);


% The genotype
genotype_list = ["WT";   %  1 = AmpKan
                 "WT";   %  2 = KanAmp(PBS)
                 "WT";   %  3 = Kan(Amp)PBS
                 "WT";   %  4 = AmpCip
                 "WT";   %  5 = CipAmp(PBS)
                 "WT";   %  6 = Cip(Amp)PBS
                 "WT";   %  7 = CipKan
                 "WT";   %  8 = KanCip(PBS)
                 "WT";   %  9 = Kan(Cip)PBS
                 "WT";   % 10 = AmpCip1000x
                 "WT";   % 11 = CipAmp1000x
                 "WT"];  % 12 = CipKan1000x
genotype = genotype_list(experiment_index);



%% Define a color palette for plots

drug_colors = [76, 176, 9;     % Amp
               7, 60, 240;     % Cip
               224, 111, 7;     % Kan
               164, 5, 227]/255;    % PBS

% Create a dict to lookup what row of the drug color matrix to use for each
% drug
drug_color_idx_dict = configureDictionary("string","double");

drug_color_idx_dict("Amp") = 1;
drug_color_idx_dict("Cip") = 2;
drug_color_idx_dict("Kan") = 3;
drug_color_idx_dict("PBS") = 4;


%% Load Dataset & Extract Data

% Load data
datafile_path = strcat("Data_and_Descriptions/",datafile);
data = readmatrix(datafile_path);

% Convert time to hours
data(:,1) = data(:,1)/60;

% Determine number of cultures
n_cults = (size(data,2)-1)/2; % Look at number of columns. 
                              % Minus 1 for the time, /2 for mean & stddev

% Determine index of first timepoint in each culture
first_time_index = -1*ones(n_cults,1);
for i = 1:n_cults
    first_time_index(i) = find(~isnan(data(:,1+i)),1,'first');
end

% Determine index of last timepoint in each culture
last_time_index = -1*ones(n_cults,1);
for i = 1:n_cults
    last_time_index(i) = find(~isnan(data(:,1+i)),1,'last');
end

% Create a cell array of (time, mean cell ct, std dev cell ct) datasets for
% each dilution condition
t_mean_stddev = {};
for i = 1:n_cults
    t_mean_stddev{i} = data(first_time_index(i):last_time_index(i),[1,1+i,1+n_cults+i]);
end

%% Plot Data (Raw)

% Figure size
fig1_pixels_width=600*(5/6);
fig1_pixels_height=500*(5/6);

% Create figure and axis
f1=figure(1);
f1.Position = [0, 0, fig1_pixels_width, fig1_pixels_height];
ax1 = gca;

hold on;
% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title("Time-Kill Assay Data");
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
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
    ax1_legends{i} = drugs(1,i);
end

legend(ax1_legends,'Location','northeast','orientation','vertical');
ax1.FontSize = 14;

% Save plot
fig1_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot1_ExperimentalData.png");
saveas(f1,fig1_plot_name);

%% Plot Data (Normalized to 1)

% Figure size
fig2_pixels_width=600*(5/6);
fig2_pixels_height=500*(5/6);

% Create figure and axis
f2=figure(2);
f2.Position = [0, 0, fig2_pixels_width, fig2_pixels_height];
ax2 = gca;

hold on;
% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title("Time-Kill Assay Data");
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
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
    ax2_legends{i} = drugs(1,i);
end

legend(ax2_legends,'Location','northeast','orientation','vertical');
ax2.FontSize = 14;

% Save plot
fig2_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot1normed_ExperimentalData.png");
saveas(f2,fig2_plot_name);



%% Estimate k_1 for the first drug from its first linear region

% Truncate data to the first linear region
n_tpts = sum(t_mean_stddev{1}(:,1) <= t_k1_drug1_end);
t_mean_stddev_k1est_drug1 = t_mean_stddev{1}(1:n_tpts,:);

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k1_drug1_fits = zeros(3,1);
% STRUCTURE OF INFO:
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_1 estimate
linear_k1drug1fit_coefs_storage = zeros(1,2);
% Structure: contains the 2 polynomial coefficients of
% the linear fit to k_1. These are stored for later plotting.

% Perform linear fit and conversion to k_1 estimate
linear_fit_coefs = polyfit(t_mean_stddev_k1est_drug1(:,1),t_mean_stddev_k1est_drug1(:,2),poly_order);
linear_k1drug1fit_coefs_storage = linear_fit_coefs;

r_sqrd_mat = corrcoef(t_mean_stddev_k1est_drug1(:,2), polyval(linear_fit_coefs,t_mean_stddev_k1est_drug1(:,1)));
k1_drug1_fits(1) = r_sqrd_mat(1,2); % R_squared

E = t_mean_stddev_k1est_drug1(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k1est_drug1(:,1));
SSE = sum(E.^2);
k1_drug1_fits(2) = sqrt(SSE/length(t_mean_stddev_k1est_drug1(:,2))); % RMSE

% Multiply the slope of the linear fit by ln(10) to estimate k_1 value
k1_drug1_fits(3) = -1*linear_fit_coefs(1)*log(10); % estimated k_1
k1_drug1_est = k1_drug1_fits(3);

%% Plot Drug 1 k_1 fit (Raw)

% Figure size
fig3_pixels_width=600*(5/6);
fig3_pixels_height=500*(5/6);

% Create figure and axis
f3=figure(3);
f3.Position = [0, 0, fig3_pixels_width, fig3_pixels_height];
ax3 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k1drug1fit_coefs_storage(1,:),0:0.01:9),'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for First Drug (%s) k_1 Estimation",drugs(1,1)));
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

legend_str = strcat(drugs(1,1),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_drug1_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{1}: %.3f h^{-1}",k1_drug1_est),"FontSize",14);
ax3.FontSize = 14;

% Save plot
fig3_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot2_k1EstimationDrug1.png");
saveas(f3,fig3_plot_name);

%% Plot Drug 1 k_1 fit (Normalized to 1)

% Figure size
fig4_pixels_width=600*(5/6);
fig4_pixels_height=500*(5/6);

% Create figure and axis
f4=figure(4);
f4.Position = [0, 0, fig4_pixels_width, fig4_pixels_height];
ax4 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k1drug1fit_coefs_storage(1,:),0:0.01:9)-t_mean_stddev{1}(1,2)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for First Drug (%s) k_1 Estimation",drugs(1,1)));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

legend_str = strcat(drugs(1,1),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_drug1_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{1}: %.3f h^{-1}",k1_drug1_est),"FontSize",14);
ax4.FontSize = 14;

% Save plot
fig4_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot2normed_k1EstimationDrug1.png");
saveas(f4,fig4_plot_name);



%% Estimate k_2 for the first drug from its second linear region

% Truncate data to the second linear region
n_tpts = sum(t_mean_stddev{1}(:,1) >= t_k2_drug1_start);
t_mean_stddev_k2est_drug1 = t_mean_stddev{1}((end-n_tpts+1):end,:);

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k2_drug1_fits = zeros(3,1);
% STRUCTURE OF INFO: 
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_2 estimate
linear_k2drug1fit_coefs_storage = zeros(1,2);
% Structure: contains the 2 polynomial coefficients of
% the linear fit to k_2. These are stored for later plotting.

% Perform linear fit and conversion to k_2 estimate
linear_fit_coefs = polyfit(t_mean_stddev_k2est_drug1(:,1),t_mean_stddev_k2est_drug1(:,2),poly_order);
linear_k2drug1fit_coefs_storage=linear_fit_coefs;

r_sqrd_mat = corrcoef(t_mean_stddev_k2est_drug1(:,2), polyval(linear_fit_coefs,t_mean_stddev_k2est_drug1(:,1)));
k2_drug1_fits(1) = r_sqrd_mat(1,2); % R_squared

E = t_mean_stddev_k2est_drug1(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k2est_drug1(:,1));
SSE = sum(E.^2);
k2_drug1_fits(2) = sqrt(SSE/length(t_mean_stddev_k2est_drug1(:,2))); % RMSE

% Multiply the slope of the linear fit by ln(10) to estimate k_2 value
k2_drug1_fits(3) = -1*linear_fit_coefs(1)*log(10); % estimated k_2
k2_drug1_est = k2_drug1_fits(3);

%% Plot Drug 1 k_2 fit (Raw)

% Figure size
fig5_pixels_width=600*(5/6);
fig5_pixels_height=500*(5/6);

% Create figure and axis
f5=figure(5);
f5.Position = [0, 0, fig5_pixels_width, fig5_pixels_height];
ax5 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k2drug1fit_coefs_storage(1,:),0:0.01:9),'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for First Drug (%s) k_2 Estimation",drugs(1,1)));
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

legend_str = strcat(drugs(1,1),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_drug1_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{2}: %.3f h^{-1}",k2_drug1_est),"FontSize",14);
ax5.FontSize = 14;

% Save plot
fig5_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot3_k2EstimationDrug1.png");
saveas(f5,fig5_plot_name);

%% Plot Drug 1 k_2 fit (Normalized to 1)

% Figure size
fig6_pixels_width=600*(5/6);
fig6_pixels_height=500*(5/6);

% Create figure and axis
f6=figure(6);
f6.Position = [0, 0, fig6_pixels_width, fig6_pixels_height];
ax6 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k2drug1fit_coefs_storage(1,:),0:0.01:9)-t_mean_stddev{1}(1,2)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for First Drug (%s) k_2 Estimation",drugs(1,1)));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

legend_str = strcat(drugs(1,1),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_drug1_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{2}: %.3f h^{-1}",k2_drug1_est),"FontSize",14);
ax6.FontSize = 14;

% Save plot
fig6_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot3normed_k2EstimationDrug1.png");
saveas(f6,fig6_plot_name);



%% Estimate k_1 for the second drug from its first linear region

% Truncate data to the first linear region
n_tpts = sum(t_mean_stddev{2}(:,1) <= t_k1_drug2_end);
t_mean_stddev_k1est_drug2 = t_mean_stddev{2}(1:n_tpts,:);

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k1_drug2_fits = zeros(3,1);
% STRUCTURE OF INFO:
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_1 estimate
linear_k1drug2fit_coefs_storage = zeros(1,2);
% Structure: contains the 2 polynomial coefficients of
% the linear fit to k_1. These are stored for later plotting.

% Perform linear fit and conversion to k_1 estimate
linear_fit_coefs = polyfit(t_mean_stddev_k1est_drug2(:,1),t_mean_stddev_k1est_drug2(:,2),poly_order);
linear_k1drug2fit_coefs_storage = linear_fit_coefs;

r_sqrd_mat = corrcoef(t_mean_stddev_k1est_drug2(:,2), polyval(linear_fit_coefs,t_mean_stddev_k1est_drug2(:,1)));
k1_drug2_fits(1) = r_sqrd_mat(1,2); % R_squared

E = t_mean_stddev_k1est_drug2(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k1est_drug2(:,1));
SSE = sum(E.^2);
k1_drug2_fits(2) = sqrt(SSE/length(t_mean_stddev_k1est_drug2(:,2))); % RMSE

% Multiply the slope of the linear fit by ln(10) to estimate k_1 value
k1_drug2_fits(3) = -1*linear_fit_coefs(1)*log(10); % estimated k_1
k1_drug2_est = k1_drug2_fits(3);

%% Plot Drug 2 k_1 fit (Raw)

% Figure size
fig7_pixels_width=600*(5/6);
fig7_pixels_height=500*(5/6);

% Create figure and axis
f7=figure(7);
f7.Position = [0, 0, fig7_pixels_width, fig7_pixels_height];
ax7 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k1drug2fit_coefs_storage(1,:),0:0.01:9),'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for Second Drug (%s) k_1 Estimation",drugs(1,2)));
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

legend_str = strcat(drugs(1,2),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_drug2_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{1}: %.3f h^{-1}",k1_drug2_est),"FontSize",14);
ax7.FontSize = 14;

% Save plot
fig7_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4_k1EstimationDrug2.png");
saveas(f7,fig7_plot_name);

%% Plot Drug 2 k_1 fit (Normalized to 1)

% Figure size
fig8_pixels_width=600*(5/6);
fig8_pixels_height=500*(5/6);

% Create figure and axis
f8=figure(8);
f8.Position = [0, 0, fig8_pixels_width, fig8_pixels_height];
ax8 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k1drug2fit_coefs_storage(1,:),0:0.01:9)-t_mean_stddev{1}(1,2)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for Second Drug (%s) k_1 Estimation",drugs(1,2)));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

legend_str = strcat(drugs(1,2),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k1_drug2_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{1}: %.3f h^{-1}",k1_drug2_est),"FontSize",14);
ax8.FontSize = 14;

% Save plot
fig8_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot4normed_k1EstimationDrug2.png");
saveas(f8,fig8_plot_name);



%% Estimate k_2 for the second drug from its second linear region

% Truncate data to the second linear region
n_tpts = sum(t_mean_stddev{2}(:,1) >= t_k2_drug2_start);
t_mean_stddev_k2est_drug2 = t_mean_stddev{2}((end-n_tpts+1):end,:);

% Linear fits to the truncated data and calculate correlation coefficient
poly_order = 1; % Linear fit is a 1st order polynomial

% Store fitting results
k2_drug2_fits = zeros(3,1);
% STRUCTURE OF INFO: 
% the first row's entry is the R_squared value,
% the second row's entry is the RMSE, 
% the third row's entry is the k_2 estimate
linear_k2drug2fit_coefs_storage = zeros(1,2);
% Structure: contains the 2 polynomial coefficients of
% the linear fit to k_2. These are stored for later plotting.

% Perform linear fit and conversion to k_2 estimate
linear_fit_coefs = polyfit(t_mean_stddev_k2est_drug2(:,1),t_mean_stddev_k2est_drug2(:,2),poly_order);
linear_k2drug2fit_coefs_storage=linear_fit_coefs;

r_sqrd_mat = corrcoef(t_mean_stddev_k2est_drug2(:,2), polyval(linear_fit_coefs,t_mean_stddev_k2est_drug2(:,1)));
k2_drug2_fits(1) = r_sqrd_mat(1,2); % R_squared

E = t_mean_stddev_k2est_drug2(:,2) - polyval(linear_fit_coefs,t_mean_stddev_k2est_drug2(:,1));
SSE = sum(E.^2);
k2_drug2_fits(2) = sqrt(SSE/length(t_mean_stddev_k2est_drug2(:,2))); % RMSE

% Multiply the slope of the linear fit by ln(10) to estimate k_2 value
k2_drug2_fits(3) = -1*linear_fit_coefs(1)*log(10); % estimated k_2
k2_drug2_est = k2_drug2_fits(3);

%% Plot Drug 2 k_2 fit (Raw)

% Figure size
fig9_pixels_width=600*(5/6);
fig9_pixels_height=500*(5/6);

% Create figure and axis
f9=figure(9);
f9.Position = [0, 0, fig9_pixels_width, fig9_pixels_height];
ax9 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k2drug2fit_coefs_storage(1,:),0:0.01:9),'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for Second Drug (%s) k_2 Estimation",drugs(1,2)));
xlabel("Time (h)");
ylabel("log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

legend_str = strcat(drugs(1,2),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_drug2_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,5.75,sprintf("Estimated k_{2}: %.3f h^{-1}",k2_drug2_est),"FontSize",14);
ax9.FontSize = 14;

% Save plot
fig9_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot5_k2EstimationDrug2.png");
saveas(f9,fig9_plot_name);

%% Plot Drug 2 k_2 fit (Normalized to 1)

% Figure size
fig10_pixels_width=600*(5/6);
fig10_pixels_height=500*(5/6);

% Create figure and axis
f10=figure(10);
f10.Position = [0, 0, fig10_pixels_width, fig10_pixels_height];
ax10 = gca;

hold on;
% Plot the linear fits as dotted lines
plot(0:0.01:9,polyval(linear_k2drug2fit_coefs_storage(1,:),0:0.01:9)-t_mean_stddev{1}(1,2)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','--','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'o','color','k');
end

title(sprintf("Fitting for Second Drug (%s) k_2 Estimation",drugs(1,2)));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end
grid on;

legend_str = strcat(drugs(1,2),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",k2_drug2_fits(2)));
legend(legend_str,'Location','northeast','orientation','vertical');
text(3.75,0,sprintf("Estimated k_{2}: %.3f h^{-1}",k2_drug2_est),"FontSize",14);
ax10.FontSize = 14;

% Save plot
fig10_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot5normed_k2EstimationDrug2.png");
saveas(f10,fig10_plot_name);


%% Fitting curves to experimental data

% Since the sequential drugging experiments do not have multiple seed
% culture dilution factors, f_1 for Drug 1 cannot be fit independently of
% P_0 for Drug 1. Similarly, there is no way to independently estimate P_0
% and f_1 for Drug 2.

% Also consider that the purpose of the sequential drugging experiments is
% to understand whether or not cells that are persistent to Drug 1 are also
% persistent to Drug 2. The purpose of the sequential drugging experiments
% is not to study the drug-induced persister formation rate.

% Therefore, the sequential drugging fits will assume f_1 = 0 for both 
% drugs. 

% P_0 for Drug 1 will be estimated by back-extrapolating the k_2
% region of the Drug 1 curve to time t=0. P_"0" for Drug 2 will be 
% estimated by back-extrapolating the k_2 region of the Drug 2 curve to the 
% time at which the second drug was added. In botyh cases, the back-
% extrapolated initial persisters are made up of an unknown ratio of 
% pre-existing and drug-induced persisters.

% Back-extrapolated initial persister CFU/mL for Drug 1:
P_0_1 = 10^polyval(linear_k2drug1fit_coefs_storage(1,:),t_mean_stddev{1}(1,1));

% Back-extrapolated initial persister CFU/mL for Drug 2:
P_0_2 = 10^polyval(linear_k2drug2fit_coefs_storage(1,:),t_mean_stddev{2}(1,1));

%% Simulate the Drug 1 curve's dynamics

% Set simulation start and end times
    t_0 = t_mean_stddev{1}(1,1); % hours
    t_end = t_mean_stddev{1}(end,1); % hours
    t_span = [t_0, t_end];

% Calculate simulation initial conditions
    N_0 = 10^t_mean_stddev{1}(1,2); % Total cells
    P_0 = P_0_1; % Back-extrapolated persister cells
    S_0 = N_0-P_0; % Susceptible cells
    init_conds = [S_0, P_0];

% Run simulation
    ode_soln_struct_Drug1 = ode45(@(t,y) model_ODEs(t,y,k1_drug1_est,k2_drug1_est,0), t_span, init_conds);

% Calculate RMSE of fit
    Drug_1_sim_vals = deval(ode_soln_struct_Drug1,t_mean_stddev{1}(:,1));
    Drug_1_sim_vals_Ntot = sum(Drug_1_sim_vals,1);
    Drug_1_sim_RMSE = rmse(log10(Drug_1_sim_vals_Ntot)',t_mean_stddev{1}(:,2));

% Extract simulation data to plot
    t_ode_Drug1 = ode_soln_struct_Drug1.x; % hours
    log10_SPN_ode_Drug1 = log10([ode_soln_struct_Drug1.y(1,:)',ode_soln_struct_Drug1.y(2,:)',sum(ode_soln_struct_Drug1.y,1)']);


%% Simulate the Drug 2 curve's dynamics

% Set simulation start and end times
    t_0 = t_mean_stddev{2}(1,1); % hours
    t_end = t_mean_stddev{2}(end,1); % hours
    t_span = [t_0, t_end];

% Calculate simulation initial conditions
    % Take the total number of cells from the Drug 1 simulation
    SP_at_tDrug2 = deval(ode_soln_struct_Drug1,t_mean_stddev{2}(1,1));
    N_0 = sum(SP_at_tDrug2);

    % Take the initial number of persisters from the back-extrapolation of
    % Drug 2's k_2 linear fit.
    P_0 = P_0_2; % Back-extrapolated persister cells

    % If the estimate for the initial number of persisters exceeds the
    % estimate for the initial total number of cells, round the persister
    % estimate down to the total cell estimate
    if P_0 > N_0
        P_0 = N_0;
    end
    
    % Calculate susceptible cells
    S_0 = N_0-P_0;

    % Set initial conditions
    init_conds = [S_0, P_0];

% Run simulation
    ode_soln_struct_Drug2 = ode45(@(t,y) model_ODEs(t,y,k1_drug2_est,k2_drug2_est,0), t_span, init_conds);

% Calculate RMSE of fit
    Drug_2_sim_vals = deval(ode_soln_struct_Drug2,t_mean_stddev{2}(:,1));
    Drug_2_sim_vals_Ntot = sum(Drug_2_sim_vals,1);
    Drug_2_sim_RMSE = rmse(log10(Drug_2_sim_vals_Ntot)',t_mean_stddev{2}(:,2));

% Extract simulation data to plot
    t_ode_Drug2 = ode_soln_struct_Drug2.x; % hours
    log10_SPN_ode_Drug2 = log10([ode_soln_struct_Drug2.y(1,:)',ode_soln_struct_Drug2.y(2,:)',sum(ode_soln_struct_Drug2.y,1)']);

%% Plot final simulation results (Raw)

% Figure size
fig11_pixels_width=600*(5/6);
fig11_pixels_height=500*(5/6);

% Create figure and axis
f11=figure(11);
f11.Position = [0, 0, fig11_pixels_width, fig11_pixels_height];
ax11 = gca;

hold on;

% Plot the ODE simulation results (N)
plot(t_ode_Drug1,log10_SPN_ode_Drug1(:,3),'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','-','LineWidth',1);
plot(t_ode_Drug2,log10_SPN_ode_Drug2(:,3),'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','-','LineWidth',1);

% Scatter plot of mean log10 CFU/mL values
point_size = 25;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),point_size,'MarkerEdgeColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',2,'MarkerEdgeAlpha',0.75);
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2),t_mean_stddev{i}(:,3),'o','color','k');
end

title(strcat("Fit to Sequential Drug Data (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Mean log_{10}(CFU/mL)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-2,9]);
    yticks(-2:9);
else
    ylim([1,9]);
    yticks(1:9);
end
grid on;

% Generate legend entries
ax11_legends = cell(n_cults,1);
ax11_legends{1} = strcat(drugs(1,1),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",Drug_1_sim_RMSE));
ax11_legends{2} = strcat(drugs(1,2),", RMSE=",sprintf("%.3f log_{10}(CFU/mL)",Drug_2_sim_RMSE));

legend(ax11_legends,'Location','northeast','orientation','vertical');

ax11.FontSize = 14;

% Save plot
fig11_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot6_FinalFit.png");
saveas(f11,fig11_plot_name);

%% Plot final f_1 and P_0 fitting results (Normalized to 1)

% Figure size
fig12_pixels_width=600*(5/6);
fig12_pixels_height=500*(5/6);

% Create figure and axis
f12=figure(12);
f12.Position = [0, 0, fig12_pixels_width, fig12_pixels_height];
ax12 = gca;

hold on;
% Plot the ODE simulation results (N)
plot(t_ode_Drug1,log10_SPN_ode_Drug1(:,3)-log10_SPN_ode_Drug1(1,3)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,1)),:),'LineStyle','-','LineWidth',3);
% t_DRUG1_sim = t_ode_Drug1;
% y_DRUG1_sim = log10_SPN_ode_Drug1(:,3)-log10_SPN_ode_Drug1(1,3)+2;
% save("t_DRUG1_sim.mat","t_DRUG1_sim");
% save("y_DRUG1_sim.mat","y_DRUG1_sim");
plot(t_ode_Drug2,log10_SPN_ode_Drug2(:,3)-log10_SPN_ode_Drug1(1,3)+2,'Color',drug_colors(drug_color_idx_dict(drugs(1,2)),:),'LineStyle','-','LineWidth',3);
% t_DRUG2_sim = t_ode_Drug2;
% y_DRUG2_sim = log10_SPN_ode_Drug2(:,3)-log10_SPN_ode_Drug1(1,3)+2;
% save("t_DRUG2_sim.mat","t_DRUG2_sim");
% save("y_DRUG2_sim.mat","y_DRUG2_sim");

% Scatter plot of mean log10 CFU/mL values
point_size = 100;
for i = 1:n_cults
    scatter(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,point_size,'filled','MarkerEdgeColor','k','MarkerFaceColor',drug_colors(drug_color_idx_dict(drugs(1,i)),:),'LineWidth',1,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
    if i == 1
        % t_DRUG1_exper = t_mean_stddev{i}(:,1);
        % y_DRUG1_exper = t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2;
        % save("t_DRUG1_exper.mat","t_DRUG1_exper");
        % save("y_DRUG1_exper.mat","y_DRUG1_exper");
    elseif i == 2
        % t_DRUG2_exper = t_mean_stddev{i}(:,1);
        % y_DRUG2_exper = t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2;
        % save("t_DRUG2_exper.mat","t_DRUG2_exper");
        % save("y_DRUG2_exper.mat","y_DRUG2_exper");
    end
end

% Put +/- 1 std dev error bars
for i = 1:n_cults
    errorbar(t_mean_stddev{i}(:,1),t_mean_stddev{i}(:,2)-t_mean_stddev{1}(1,2)+2,t_mean_stddev{i}(:,3),'color','k','LineWidth',1,'LineStyle','none','CapSize',12);
    if i == 1
        % SD_DRUG1_exper = t_mean_stddev{i}(:,3);
        % save("SD_DRUG1_exper.mat","SD_DRUG1_exper");
    elseif i == 2
        % SD_DRUG2_exper = t_mean_stddev{i}(:,3);
        % save("SD_DRUG2_exper.mat","SD_DRUG2_exper");
    end
end

%title(strcat("Fit to Sequential Drug Data (",experiment_nickname,")"));
xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
if (drugs(1)=="Kan" || drugs(2)=="Kan")
    ylim([-8,2.2]);
    yticks(-8:1:2);
    yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
else
    ylim([-5,2.2]);
    yticks(-5:1:2);
    yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
end

grid on;
set(ax12,'GridAlpha',0.4);

% Generate legend entries
ax12_legends = cell(n_cults,1);
ax12_legends{1} = drugs(1,1);
ax12_legends{2} = drugs(1,2);

legend(ax12_legends,'Location','northeast','orientation','vertical');

ax12.FontSize = 28;

% Save plot
fig12_plot_name = strcat("Output_Plots/",experiment_nickname,"_Plot6normed_FinalFit.png");
saveas(f12,fig12_plot_name);




end




%% Functions

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end