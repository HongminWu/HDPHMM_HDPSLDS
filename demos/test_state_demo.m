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
 for testing_trial =2:21 % for testing each state
    close all;
    %channel-#1
    static_test_state(testing_trial);  % push the whole data as the observation, and get the results directly.
    
    %channel-#2
    %dynamic_test_state(testing_trial); %Obtained the results by dynamical display the observation
 end
 cd (rootPath)
toc;
end