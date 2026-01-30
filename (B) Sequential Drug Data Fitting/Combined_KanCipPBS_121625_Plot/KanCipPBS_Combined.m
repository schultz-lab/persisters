%% Clear
clear;
clc;
close all;


%% Define a color palette for plots

drug_colors = [76, 176, 9;     % Amp
               7, 60, 240;     % Cip
               224, 111, 7;     % Kan
               164, 5, 227]/255;    % PBS

%% Load the data to plot (saved from the main sequential fitting & plotting code)

load("t_Kan_sim.mat");
load("y_Kan_sim.mat");

load("t_PBS_sim.mat");
load("y_PBS_sim.mat");

load("t_Cip_sim.mat");
load("y_Cip_sim.mat");


load("t_Kan_exper.mat");
load("y_Kan_exper.mat");
load("SD_Kan_exper.mat");

load("t_PBS_exper.mat");
load("y_PBS_exper.mat");
load("SD_PBS_exper.mat");

load("t_Cip_exper.mat");
load("y_Cip_exper.mat");
load("SD_Cip_exper.mat");


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
plot(t_Kan_sim,y_Kan_sim,'Color',drug_colors(3,:),'LineStyle','-','LineWidth',3);
plot(t_PBS_sim,y_PBS_sim,'Color',drug_colors(4,:),'LineStyle','-','LineWidth',3);
plot(t_Cip_sim,y_Cip_sim,'Color',drug_colors(2,:),'LineStyle','-','LineWidth',3);

% Scatter plot of mean log10 CFU/mL values
point_size = 100;
scatter(t_Kan_exper,y_Kan_exper,point_size,'filled','MarkerEdgeColor','k','MarkerFaceColor',drug_colors(3,:),'LineWidth',1,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
scatter(t_PBS_exper,y_PBS_exper,point_size,'filled','MarkerEdgeColor','k','MarkerFaceColor',drug_colors(4,:),'LineWidth',1,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);
scatter(t_Cip_exper,y_Cip_exper,point_size,'filled','MarkerEdgeColor','k','MarkerFaceColor',drug_colors(2,:),'LineWidth',1,'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.75);

% Put +/- 1 std dev error bars
errorbar(t_Kan_exper,y_Kan_exper,SD_Kan_exper,'color','k','LineWidth',1,'LineStyle','none','CapSize',12);
errorbar(t_PBS_exper,y_PBS_exper,SD_PBS_exper,'color','k','LineWidth',1,'LineStyle','none','CapSize',12);
errorbar(t_Cip_exper,y_Cip_exper,SD_Cip_exper,'color','k','LineWidth',1,'LineStyle','none','CapSize',12);


xlabel("Time (h)");
ylabel("Percent Survival (%)");
xlim([0,8]);
xticks(0:8);
ylim([-8,2.2]);
yticks(-8:1:2);
yticklabels(["10^{-8}","10^{-7}","10^{-6}","10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
grid on;
set(ax12,'GridAlpha',0.4);

% Generate legend entries
ax12_legends = cell(3,1);
ax12_legends{1} = "Kan";
ax12_legends{2} = "PBS";
ax12_legends{3} = "Cip";

legend(ax12_legends,'Location','northeast','orientation','vertical');

ax12.FontSize = 28;

% Save plot
saveas(f12,"KanCipPBS_Combined.png");