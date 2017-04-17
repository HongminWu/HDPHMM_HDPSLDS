function demo_testing()
clc;
clear;
clear global;
    global_variables;
    disp(strcat('Learned_Models:', modelPath));
    disp('.');
    disp('........................................');
    disp('.');
    disp(strcat('Testing_Observations:', TESTING_DATASET_PATH));
 for testing_trial = 7:7
    close all;
    %testing_trial = 1;
    testing(testing_trial)
 end
end