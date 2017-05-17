clear global;
clc; 
global rootPath rootDATApath modelPath                            %path
global METHOD ROBOT TASK STATE SIGNAL_TYPE TRUNCATION_STATES TRUNCATION_COMPONENTS %experiment set up
global TRAINING_SIM_REAL TRAINING_SUCCESS_FAILURE TRAINING_DATASET TRAINING_DATASET_PATH TRAINING_PLOT_ON PLOT_SAVE saveDir Niter %training
global TESTING_SIM_REAL TESTING_SUCCESS_FAILURE TESTING_DATASET TESTING_PLOT_ON TESTING_DATASET_PATH TESTING_RESULTS_PATH %testing
global THRESHOLD_PATH DIS_PERIOD TIME_STEP WRENCH_DERIVATIVE
global TRUE_POSITIVE  FALSE_NEGATIVE FALSE_POSITIVE TRUE_NEGATIVE
% general initialization
TRUNCATION_STATES           = 20;    %states : truncation level of the DP(Dirichlet Process) prior on HMM transition distributions pi_k
TRUNCATION_COMPONENTS       = 1;     %components of mixture gaussian: truncation level of the DPMM(Dirichlet Process Mixture Model) on emission distributions pi_s
DIS_PERIOD                  = 10;    %define the display period for echoing the the screen
TIME_STEP                   = 0.005; %time-step for recording the signals
METHOD                      = 'HDPVARHMM(2)'; %sHDPHMM, HDPVARHMM(1),HDPVARHMM(2)
ROBOT                       = 'HIRO';
TASK                        = 'SA';
STATE                       = {'APPROACH', 'ROTATION', 'INSERTION', 'MATING'};
SIGNAL_TYPE                 = {'R_Torques'}; %'R_CartPos'  %'R_CartPos_Corrected' %'R_Torques'  %'R_Angles'
WRENCH_DERIVATIVE           = true;
rootPath                    = '/home/Homls/npBayesHMM/matlab/HDPHMM_HDPSLDS';
rootDATApath                = '/home/Homls/npBayesHMM/HIRO_SA_DATA';
multiModal                  = ''; 
for nSignal                 = 1 : length(SIGNAL_TYPE) 
    multiModal              = strcat(multiModal,'_',SIGNAL_TYPE{nSignal}); 
end

addpath(genpath(rootPath));
addpath(genpath(rootDATApath));

% only for training
TRAINING_SIM_REAL           = 'REAL'; %REAL SIM
TRAINING_SUCCESS_FAILURE    = 'SUCCESS'; %SUCCESS
TRAINING_DATASET            = 'REAL_HIRO_ONE_SA_SUCCESS';%REAL_HIRO_ONE_SA_SUCCESS SIM_HIRO_ONE_SA_SUCCESS
TRAINING_DATASET_PATH       = strcat(rootDATApath,'/',TRAINING_DATASET);
modelPath                   = strcat(rootPath,'/learned_models','/',METHOD,'_',TRAINING_SIM_REAL,'_',TRAINING_SUCCESS_FAILURE,'_',ROBOT,'_',TASK,multiModal,'_');
saveDir                     = strcat(rootPath,'/saveDir');
TRAINING_PLOT_ON            = false; %true
Niter                       = 500;

% only for testing
TESTING_SIM_REAL            = 'REAL'; %REAL SIM
TESTING_SUCCESS_FAILURE     = 'FAILURE'; % SUCCESS FAILURE
TESTING_DATASET             = 'REAL_HIRO_ONE_SA_ERROR_CHARAC'; %REAL_HIRO_ONE_SA_SUCCESS, REAL_HIRO_ONE_SA_ERROR_CHARAC
TESTING_DATASET_PATH        = strcat(rootDATApath,'/',TESTING_DATASET);
TESTING_RESULTS_PATH        = strcat(rootPath,'/testing_results','/',METHOD,'_',TESTING_SIM_REAL,'_',TESTING_SUCCESS_FAILURE,multiModal);
TESTING_PLOT_ON             = true;
PLOT_SAVE                   = true;

%only for calculate state_threshold
THRESHOLD_PATH              = strcat(rootPath,'/state_threshold/',METHOD,'_',multiModal,'_');

% for ROC curve
TRUE_POSITIVE     = 0;  % SUCCESS -> predicted as -> SUCCESS
FALSE_NEGATIVE = 0;  % SUCCESS -> predicted as ->  FAILURE
FALSE_POSITIVE  = 0;  %  FAILURE -> predicted as ->  SUCCESS 
TRUE_NEGATIVE  = 0;  %  FAILURE -> predicted as ->  FAILURE 


