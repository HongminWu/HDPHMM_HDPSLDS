% for k-fold cross validation
function cross_validation()
clc;
clear;
clear global;
global TESTING_SIM_REAL TESTING_SUCCESS_FAILURE TESTING_DATASET 
global rootPath  rootDATApath TESTING_DATASET_PATH TESTING_RESULTS_PATH   
global TRUE_POSITIVE  FALSE_NEGATIVE FALSE_POSITIVE TRUE_NEGATIVE


succIdx1  = {'02','03','04','05','06','07','08','09','10','11'};
succIdx2  = {'12','13','14','15','16','17','18','19','20','21'};
failIdx1  = {'02','03','04','05','06','07','08','09'};
failIdx2  = {'10','11','12','13','14','15','16','17'};
succ = {succIdx1,succIdx2};
fail = {failIdx1,failIdx2};
%2-fold
%     comb_1    ={succIdx1,failIdx1};
%     comb_2    ={succIdx1, failIdx2};
%     comb_3    ={succIdx2,failIdx1};
%     comb_4    ={succIdx2,failIdx2};

% succIdx1  = {'02','03','04','05','06'};
% succIdx2  =  {'07','08','09','10','11'};
% succIdx3  = {'12','13','14','15','16'};
% succIdx4  = {'17','18','19','20','21'};
% failIdx1  = {'02','03','04','05'};
% failIdx2=  {'06','07','08','09'};
% failIdx3  = {'10','11','12','13'};
% failIdx4 = {'14','15','16','17'};
% succ = {succIdx1,succIdx2,succIdx3,succIdx4};
% fail =   {failIdx1,failIdx2,failIdx3,failIdx4};

nCross   = 0;
ROC      = [];
for idx = 1:length(succ)     
    calculate_state_threshold(succ{idx}); % Calculate the threshold and generate the optimal model among the training       
    global_variables;
    for jdx = 1:length(fail)
        nCross = nCross +1; 
        for sIdx = 1:length(succ{jdx})        %testing success trials
            % only for testing
            TESTING_SIM_REAL            = 'REAL'; %REAL SIM
            TESTING_SUCCESS_FAILURE     = 'SUCCESS'; % SUCCESS FAILURE
            TESTING_DATASET             = 'REAL_HIRO_ONE_SA_SUCCESS'; %REAL_HIRO_ONE_SA_SUCCESS, REAL_HIRO_ONE_SA_ERROR_CHARAC
            TESTING_DATASET_PATH        = strcat(rootDATApath,'/',TESTING_DATASET);
            TESTING_RESULTS_PATH        = strcat(rootPath,'/Cross_Validtation','/',num2str(nCross));
            static_test_state(char(succ{jdx}(sIdx)));
        end
    
        close all;
        
        for fIdx = 1:length(fail{jdx})       %testing success trials
            % only for testing
            TESTING_SIM_REAL            = 'REAL'; %REAL SIM
            TESTING_SUCCESS_FAILURE     = 'FAILURE'; % SUCCESS FAILURE
            TESTING_DATASET             = 'REAL_HIRO_ONE_SA_ERROR_CHARAC'; %REAL_HIRO_ONE_SA_SUCCESS, REAL_HIRO_ONE_SA_ERROR_CHARAC
            TESTING_DATASET_PATH        = strcat(rootDATApath,'/',TESTING_DATASET);
            TESTING_RESULTS_PATH        = strcat(rootPath,'/Cross_Validtation','/',num2str(nCross));
            static_test_state(char(fail{jdx}(fIdx)));
        end

        TPR =  (TRUE_POSITIVE/(TRUE_POSITIVE + FALSE_NEGATIVE));
        FPR = (FALSE_POSITIVE/(FALSE_POSITIVE + TRUE_NEGATIVE));
        ROC =[ROC; [FPR, TPR ]]
        TRUE_POSITIVE     = 0;  % SUCCESS -> predicted as -> SUCCESS
        FALSE_NEGATIVE  = 0;  % SUCCESS -> predicted as ->  FAILURE
        FALSE_POSITIVE   = 0;  %  FAILURE -> predicted as ->  SUCCESS 
        TRUE_NEGATIVE   = 0;  %  FAILURE -> predicted as ->  FAILURE
    end
end

figure;
ROC = sortrows(ROC, 1);
plot(ROC(:,1), ROC(:,2), 'r-','LineWidth', 3);
title('ROC by HDP-VAR-HMM');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
end
