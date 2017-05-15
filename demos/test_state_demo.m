function test_state_demo()
clc;
clear;
clear global;
tic;
    global_variables;
    disp(strcat('Learned_Models:', modelPath));
    disp('.');
    disp('........................................');
    disp('.');
    disp(strcat('Testing_Observations:', TESTING_DATASET_PATH));
 for testing_trial =7:7 % for testing each state
    close all;
    test_state(testing_trial)
 end
 cd (rootPath)
toc;
end