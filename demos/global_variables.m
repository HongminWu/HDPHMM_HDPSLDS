clear;
clc; 
global rootPath rootDATApath modelPath                            %path
global METHOD  ROBOT TASK STATE SIGNAL_TYPE TRUNCATION_STATES TRUNCATION_COMPONENTS %experiment set up
global TRAINING_SIM_REAL TRAINING_SUCCESS_FAILURE TRAINING_DATASET TRAINING_DATASET_PATH TRAINING_PLOT_ON saveDir Niter %training
global TESTING_SIM_REAL TESTING_SUCCESS_FAILURE TESTING_DATASET TESTING_PLOT_ON TESTING_DATASET_PATH TESTING_RESULTS_PATH LoglikelihoodRange COLOR %testing

% general initialization
TRUNCATION_STATES           = 20; %states : truncation level of the DP(Dirichlet Process) prior on HMM transition distributions pi_k
TRUNCATION_COMPONENTS       = 1;  %components of mixture gaussian: truncation level of the DPMM(Dirichlet Process Mixture Model) on emission distributions pi_s
METHOD                      = 'HDPVARHMM'; %sHDPHMM, HDPVARHMM
ROBOT                       = 'HIRO';
TASK                        = 'SA';
STATE                       = {'APPROACH', 'ROTATION', 'INSERTION', 'MATING'};
SIGNAL_TYPE                 = {'R_Torques'}; %'R_CartPos'  %'R_CartPos_Corrected' %'R_Torques'  %'R_Angles'
rootPath                    = 'F:\Matlab\Emily_Fox\HDPHMM_HDPSLDS';
rootDATApath                = 'F:\Matlab\Emily_Fox\HIRO_SA_DATA';


addpath(genpath(rootPath));
addpath(genpath(rootDATApath));

% only for training
TRAINING_SIM_REAL           = 'REAL';
TRAINING_SUCCESS_FAILURE    = 'SUCCESS'; %for saving
TRAINING_DATASET            = 'REAL_HIRO_ONE_SA_SUCCESS';%REAL_HIRO_ONE_SA_SUCCESS
TRAINING_DATASET_PATH       = strcat(rootDATApath,'/',TRAINING_DATASET);
modelPath                   = strcat(rootPath,'/learned_models','/',METHOD,'_',TRAINING_SIM_REAL,'_',TRAINING_SUCCESS_FAILURE,'_',ROBOT,'_',TASK,'_');
saveDir                     = strcat(rootPath,'/saveDir');
TRAINING_PLOT_ON            = true; %true
Niter                       = 500;

% only for testing
TESTING_SIM_REAL            = 'REAL';
TESTING_SUCCESS_FAILURE     = 'SUCCESS'; % SUCCESS FAILURE
TESTING_DATASET             = 'REAL_HIRO_ONE_SA_SUCCESS'; %REAL_HIRO_ONE_SA_SUCCESS, REAL_HIRO_ONE_SA_ERROR_CHARAC
TESTING_DATASET_PATH        = strcat(rootDATApath,'/',TESTING_DATASET);
TESTING_RESULTS_PATH        = strcat(rootPath,'/testing_results');
COLOR                       = ['r', 'g', 'b','k'];
TESTING_PLOT_ON             = true;
LoglikelihoodRange          = [-1.0000e+6, 0];  %for plot
