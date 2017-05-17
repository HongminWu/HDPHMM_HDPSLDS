%use for plot REAL_HIRO_ONE_SA_ERROR_CHARAC 
function plot_dataset()
clc;    
clear;
close all;
clear global;
global_variables;
global rootDATApath
dataPath = strcat(rootDATApath,'/','REAL_HIRO_ONE_SA_ERROR_CHARAC');
trialIdx = {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17'};
for idx = 1:length(trialIdx)
    [data, R_State, foldname] = load_one_trial(trialIdx{idx}, dataPath);
    plot(R_State');
    hold on;
end
end