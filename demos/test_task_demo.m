clc;
clear;
clear global;
    global_variables;
    disp(strcat('Learned_Models:', modelPath));
    disp('.');
    disp('........................................');
    disp('.');
    disp(strcat('Testing_Observations:', TESTING_DATASET_PATH));
    %trialIdx = {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
    trialIdx = {'06'};
     for  idx=1:length(trialIdx) % for testing each state
        close all;
        test_task(trialIdx{idx})
     end
 cd (rootPath)