function test_task_demo()
clc;
clear;
clear global;
    global_variables;
    disp(strcat('Learned_Models:', modelPath));
    disp('.');
    disp('........................................');
    disp('.');
    disp(strcat('Testing_Observations:', TESTING_DATASET_PATH));
 for testing_trial = 3:8 % for testing each state
    close all;
    test_task(testing_trial)
 end
 cd (rootPath)
end