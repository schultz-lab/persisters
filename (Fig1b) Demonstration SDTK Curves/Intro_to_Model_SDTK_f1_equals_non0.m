%% Notes
% Whenever S (susceptible cells) and P (persister cells) are together in
% the same matrix, we will order them so S is #1 and P is #2.

%% Housekeeping
clear;

%% Set up ODE solver

t_0 = 0; % Start time of the simulation, unit: hours.
t_end = 6; % End time of the simulation, unit: hours.
t_span = [t_0, t_end];

S_0 = 10^8; % Initial susceptible cell population, unit: number of cells.

%%% Looped Parameter %%%
P_0 = [10^5, 10^4, 10^3]; % Initial persister cell population, unit: number of cells.
loop_dim_1 = length(P_0);
%%% Looped Parameter %%%

k_1 = 7; % Death rate of susceptible cells, unit: 1/h.
k_2 = 0.7; % Death rate of persister cells, unit: 1/h

% Drug-induced model
f_1 = 0.5*10^-3; % Rate of susceptible to persister induction, unit: 1/h.

%% Prepare storage of results
t_array = cell(loop_dim_1,1);
y_array = cell(loop_dim_1,1);

%% Loop ODE Solver Runs & Results Storage

for i = 1:loop_dim_1

    init_conds = [S_0, P_0(i)];

    [t_ode,y_ode] = ode45(@(t,y) model_ODEs(t,y,k_1,k_2,f_1), t_span, init_conds);

    t_array{i} = t_ode;
    y_array{i} = y_ode;

end

%% Plot Results

% Size of figure window
pixels_width=500;
pixels_height=400;

% Create figure window
f1=figure(1);
f1.Position = [0, 0, pixels_width, pixels_height];

% Define a color palette for plots
colors_RGB = [1.0, 0.0, 0.0; % red
              0.3, 0.7, 0.3; % green
              0.0, 0.0, 1.0]; % blue

% Plot Data - N
semilogy(t_array{1},100*(y_array{1}(:,1)+y_array{1}(:,2))/(y_array{1}(1,1)+y_array{1}(1,2)),'-','Color',colors_RGB(1,:),'LineWidth',3);
hold on;
semilogy(t_array{2},100*(y_array{2}(:,1)+y_array{2}(:,2))/(y_array{2}(1,1)+y_array{2}(1,2)),"Color",colors_RGB(2,:),'LineStyle','-','LineWidth',3);
semilogy(t_array{3},100*(y_array{3}(:,1)+y_array{3}(:,2))/(y_array{3}(1,1)+y_array{3}(1,2)),'-','Color',colors_RGB(3,:),'LineWidth',3);

% Plot Data - P
semilogy(t_array{1},100*(y_array{1}(:,2))/(y_array{1}(1,1)+y_array{1}(1,2)),'--','Color',colors_RGB(1,:),'LineWidth',3);
semilogy(t_array{2},100*(y_array{2}(:,2))/(y_array{2}(1,1)+y_array{2}(1,2)),'Color',colors_RGB(2,:),'LineStyle','--','LineWidth',3);
semilogy(t_array{3},100*(y_array{3}(:,2))/(y_array{3}(1,1)+y_array{3}(1,2)),'--','Color',colors_RGB(3,:),'LineWidth',3);

% Figure Labels
text(0.05, 10^2.05, "N_0","Color",'k','FontSize',24,'FontWeight','bold');
text(0.05, 10^-0.6, "P_0","Color",colors_RGB(1,:),'FontSize',22,'FontWeight','bold');
text(0.05, 10^-1.425, "P_0","Color",colors_RGB(2,:),'FontSize',22,'FontWeight','bold');
text(0.05, 10^-3.2, "P_0","Color",colors_RGB(3,:),'FontSize',22,'FontWeight','bold');
text(0.7, 10^0.4, "k_1 = 7 h^{-1}","Color",'k','FontSize',24,'FontWeight','bold');
hndl = text(3.5, 10^-1.4, "k_2 = 0.7 h^{-1}","Color",'k','FontSize',24,'FontWeight','bold');
set(hndl,'Rotation',-13);
text(0.25, 10^-4.3, "f_1 = 0.0005 h^{-1}","Color",'k','FontSize',24,'FontWeight','bold');

xlabel("Time (h)");
ylabel("Percent Survival (%)");
grid on;
set(get(f1,'CurrentAxes'),'GridAlpha',0.4,'MinorGridAlpha',0.4);
ylim([10^-5 10^2.5]);
yticks(10.^(-5:1:2));
yticklabels(["10^{-5}","10^{-4}","10^{-3}","0.01","0.1","1","10","100"]);
xlim([0,6]);
xticks(0:1:6);
% title({"Persistence Induction Model","Seed Culture Dilution"});

% Font Size
axis1 = gca;
axis1.FontSize=28;

% Save Figure
saveas(f1,"Intro_to_Model_SDTK_f1_equals_non0.png");

%% model_ODEs function

function dydt = model_ODEs(~,y,k_1,k_2,f_1) % The t input isn't used, so replace it with a "~"
    % y(1) is S
    % y(2) is P 
    dydt(1,1) = -k_1*y(1) - f_1*y(1); % dydt(1,1) is dS/dt
    dydt(2,1) = f_1*y(1) - k_2*y(2); % dydt(2,1) is dP/dt
end
