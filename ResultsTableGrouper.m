%% Note: Once you've run all of the fits in SDTK_Fitting.m, run this code
% to compile all fitting results into a single table.

clear;
clc;

for i = 1:19
    if i == 1
        table_holder = load(strcat("Output_Tables/Experiment",num2str(i),"_ResultsTable.mat"));
        supertable = table_holder.fit_results_tab;
    else
        table_holder = load(strcat("Output_Tables/Experiment",num2str(i),"_ResultsTable.mat"));
        supertable = [supertable; table_holder.fit_results_tab];
    end
end


save("Output_Tables/FittingResultsSuperTable.mat","supertable");