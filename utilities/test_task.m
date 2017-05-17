%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function test_task(trialID)
global rootPath modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH 
global TRUNCATION_STATES TRUNCATION_COMPONENTS
global TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global TESTING_PLOT_ON
%% Training Model
cd (strcat(modelPath,'TASK'));
    file                     = dir('*.mat');
    learnedModel             = load(file(1).name);                                 %for initializating the 'meanmodel'
    [data, R_State, foldname]= load_one_trial(trialID,TESTING_DATASET_PATH);
cd(rootPath);
%% testing
    gHandle_testing         = figure;
    title(strcat('Testing: The likelihood of testing trial w.r.t the learned model'));
%% Inference using training model and testing trial
    settings.Kz             = TRUNCATION_STATES;  % nStats
    settings.Ks             = TRUNCATION_COMPONENTS;  % nComponents of mixture gaussian
    data_struct             = struct;
    left                    = 0.1; 
    bottow                  = 0.1; 
    width                   = 0.8; 
    heigh                   = 0.4;
    data                    = data_preprocessing(data);  %data_preprocessing !!!!!!!!!!!!!!!!!!!!!!
    sensor                  = [];
    log_likelihood          = [];
     cd (TESTING_RESULTS_PATH);
     exlikhood = load('log_likelihood.mat');
 
obsModel            = learnedModel.obsModel;
obsModelType        = learnedModel.obsModel.type;
r = 1;
if strcmp(obsModelType,'AR')
    r = r + obsModel.r;
end
    
    
for obsIdx = r:length(data)   
    disp(strcat('Testing: obsIdx = ',num2str(obsIdx),'/',num2str(length(data))));
    
    if strcmp(obsModelType,'AR') & (obsIdx == r)
        sensor                  =  data(:, obsIdx - obsModel.r: obsIdx);
    end
    sensor                  = [sensor, data(:,obsIdx)] ;
    data_struct.obs         = sensor;
    data_struct.true_labels = ones(1,size(sensor,2));
    data_struct.test_cases  = 1;
    data_struct.blockSize   = ones(1,size(sensor,2));
    data_struct.blockEnd    = 1:size(sensor,2);
    if strcmp(obsModelType,'AR')
        [X,valid]               = makeDesignMatrix(data_struct.obs,obsModel.r);
        data_struct.obs         = data_struct.obs(:,find(valid));
        data_struct.X           = X(:,find(valid));
        data_struct.blockSize   = ones(1,size(data_struct.obs,2));
        data_struct.blockEnd    = cumsum(data_struct.blockSize);
        data_struct.true_labels = data_struct.true_labels(find(valid));
    end
    %testing the observed data
    dist_struct             = learnedModel.dist_struct; 
    theta                   = learnedModel.theta;
    obsModel                = learnedModel.obsModel;
    obsModelType            = learnedModel.obsModel.type;
    [total_log_likelihood,... 
        testing_neglog_c]   = observation_likelihood(data_struct,obsModelType,dist_struct,theta); 
    
    log_likelihood          = [log_likelihood; total_log_likelihood]; %
    
    if TESTING_PLOT_ON  & (rem(obsIdx, 50) == 0)     
        cla(gHandle_testing);
        %info: subplot('position',[left bottom width height]) range from 0.0 to 1.0
        subplot_1 = subplot(2,1,1,'position', [left bottow + heigh  width heigh],'Parent',gHandle_testing);
        %subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
        title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...
              strcat('TrainingModels: ', TRAINING_SIM_REAL,'\_', TRAINING_SUCCESS_FAILURE); ...
              strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',foldname)});
        plot(log_likelihood, 'b-','LineWidth',3,'Parent',subplot_1);
        grid on
        grid MINOR
        axis auto
        xlim([0  length(data)]);
       % ylim([LoglikelihoodRange(1), LoglikelihoodRange(2)]);  
        gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
         hold on;
         plot(exlikhood.log_likelihood, 'r--','LineWidth',3,'Parent',subplot_1);
        
        subplot_2 = subplot(2,1,2,'position', [left bottow width heigh],'Parent',gHandle_testing);
        plot(data(:,1:obsIdx)','Parent',subplot_2); % for animation
        grid on
        grid MINOR
        gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
        xlim([0  length(data)]);
        pause(0.000001);
    end
end
if (exist(TESTING_RESULTS_PATH,'dir') == 0)
    mkdir(TESTING_RESULTS_PATH);
end
cd (TESTING_RESULTS_PATH);
save('log_likelihood.mat','log_likelihood')
saveas(gHandle_testing, strcat(TESTING_SIM_REAL,'_',TESTING_SUCCESS_FAILURE,'_',foldname),'jpg');
cd (rootPath);