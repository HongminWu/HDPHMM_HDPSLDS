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
    
    %for successful dataset
    trialIdx = {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
    
    %for failure dataset 'REAL_HIRO_ONE_SA_ERROR_CHARAC'
    %Approach:  trialIdx = {'14','15','16','17'};
    %Rotation:  trialIdx = {'07','08','13','15'};
    %Insertion: trialIdx = {'09','11','17'};
    %Mating:    trialIdx = {'06','10','12','16'};
    
 for  idx=1:length(trialIdx) % for testing each state
    close all;
    %channel-#1
    static_test_state(trialIdx{idx});  % push the whole data as the observation, and get the results directly.
    
    %channel-#2
    %dynamic_test_state(testing_trial); %Obtained the results by dynamical display the observation
 end
 cd (rootPath)
toc;
end