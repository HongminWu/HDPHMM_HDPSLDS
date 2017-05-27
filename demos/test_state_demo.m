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
    %trialIdx =
    trialIdx =  {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};%,'22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45'};
    %trialIdx =  {'03'};
 for  idx=1:length(trialIdx) % for testing each state
    close all;
    %channel-#1
     static_test_state_classification(trialIdx{idx});  %get the results directly.
    %static_test_failure_detection(trialIdx{idx});  
    
    %channel-#2
    %dynamic_test_state(testing_trial); %Obtained the results by dynamical display the observation
 end
cd (rootPath)
toc;