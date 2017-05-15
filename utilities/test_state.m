%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function test_state(trialID)
global rootPath modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH 
global METHOD ROBOT TASK STATE TRUNCATION_STATES TRUNCATION_COMPONENTS
global TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global COLOR TESTING_PLOT_ON THRESHOLD_PATH DIS_PERIOD
                                          
%% Training Model
learnedModel = {};
model = [];

for nlearnedModel = 1 : length(STATE)
        %find out the optimal for each state, this step is achieved by function "calculate_state_threshold"
        thresholdPath  = strcat(THRESHOLD_PATH, char(STATE(nState)));
        if (exist(thresholdPath,'dir') == 0)
            disp('Please define the optimal model for each state before testing? Training -> Model Selection -> Testing');
            break;
        end
        cd(thresholdPath);
        model_select    = load(strcat('MODEL_SELECTION_METHOD_',char(STATE(nState)),'.mat'),'MODEL_SELECTION_METHOD');
        [~, optimalIdx] = max(model_select);
        
        % Load Model for selected state. We will compare other observations with this.
        cd (strcat(modelPath,char(STATE(nlearnedModel))));
        file         = dir('*.mat');
        meanmodel = [];    
        trained_model_name = file(optimalIdx).name;
        meanmodel    = load(trained_model_name);                                 %for initializating the 'meanmodel'
  %{      
    for nfile        = 1 : size(file,1)
        temp         = load(file(nfile).name);                             % load the parameters of all learned HMM models
        model        = [model, temp];
    end  
        invSigma     = zeros(size(model(1).theta.invSigma));
        mu           = zeros(size(model(1).theta.mu));
        pi_z         = zeros(size(model(1).dist_struct.pi_z));
        pi_init      = zeros(size(model(1).dist_struct.pi_init));
        pi_s         = zeros(size(model(1).dist_struct.pi_s));
        beta_vec     = zeros(size(model(1).dist_struct.beta_vec));
    for nmodel       = 1 : size(model,2)                                   %for calculating the mean of models
        invSigma     = invSigma + model(nmodel).theta.invSigma;
        mu           = mu       + model(nmodel).theta.mu;
        pi_z         = pi_z     + model(nmodel).dist_struct.pi_z;
        pi_init      = pi_init  + model(nmodel).dist_struct.pi_init;
        pi_s         = pi_s     + model(nmodel).dist_struct.pi_s;
        beta_vec     = beta_vec + model(nmodel).dist_struct.beta_vec;
    end
        meanmodel.theta.invSigma       = invSigma/size(model,2);
        meanmodel.theta.mu             = mu/size(model,2);
        meanmodel.dist_struct.pi_z     = pi_z/size(model,2);
        meanmodel.dist_struct.pi_init  = pi_init/size(model,2);
        meanmodel.dist_struct.pi_s     = pi_s/size(model,2);
        meanmodel.dist_struct.beta_vec = beta_vec/size(model,2);
  %}
        learnedModel                   = [learnedModel, {meanmodel}];
end

[data, R_State, foldname] = load_one_trial(trialID,TESTING_DATASET_PATH);
cd(rootPath);
%% testing
gHandle_testing = figure;
title(strcat('Testing: The likelihood of testing trial w.r.t the learned model'));
sensor = [];
%threshold = []; 
temp_log_likelihood = zeros(length(STATE),1);
total_log_likelihood = [];
%% Inference using training model and testing trial
settings.Kz             =  TRUNCATION_STATES;      % nStats
settings.Ks             =  TRUNCATION_COMPONENTS;  % nComponents of mixture gaussian
data_struct             = struct;
left                    = 0.1; 
bottow                  = 0.1; 
width                   = 0.8; 
heigh                   = 0.4;
trans                   = 0;
data                    = data_preprocessing(data); % mean = 0 and covariance = Unit matrix

obsModel            = learnedModel{nlearnedModel}.obsModel;
obsModelType        = learnedModel{nlearnedModel}.obsModel.type;
r = 1;
if strcmp(obsModelType,'AR')
    r = r + obsModel.r;
end

for obsIdx = r:length(data)     
    if rem(obsIdx, DIS_PERIOD) == 0     
        disp(strcat('Testing: obsIdx = ',num2str(obsIdx),'/',num2str(length(data))));
    end
    %load different learned model at different state
    if ismember(obsIdx - obsModel.r, R_State(:,1:end-1))
        sensor    = [];
        %threshold = [];
        gridxy(R_State(:,2:end-1),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
        trans     = trans + 1;
        if strcmp(obsModelType,'AR')
            sensor =  data(:, (obsIdx - obsModel.r): obsIdx);
        end 
        %load the threshold 
        %threshold  = load(strcat(THRESHOLD_PATH, METHOD,'_', char(STATE(trans)),'/','THRESHOLD_',char(STATE(trans)),'.mat'));
    end
    sensor                  = [sensor, data(:,obsIdx)];
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
    %testing the observed data for each learned model
    for nlearnedModel       = 1 : length(STATE)
        dist_struct         = learnedModel{nlearnedModel}.dist_struct; 
        theta               = learnedModel{nlearnedModel}.theta;
        obsModel            = learnedModel{nlearnedModel}.obsModel;
        obsModelType        = learnedModel{nlearnedModel}.obsModel.type;
        [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);       
        temp_log_likelihood(nlearnedModel,1) = testing_total_log_likelihood ; %
    end
    total_log_likelihood = [total_log_likelihood, temp_log_likelihood]; %store the results
    
    if TESTING_PLOT_ON & ((rem(obsIdx, vis_period) == 0) | ismember(obsIdx, R_State(:,end)))         
        clf(gHandle_testing);
%         subplot_1 = subplot(2,1,1,'position', [left bottow + heigh  width heigh],'Parent',gHandle_testing);
        subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
%        plot(threshold.threshold','c--');              
        for sIdx = 1 : size(total_log_likelihood,1)
            title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...                                                                                                                        
                  strcat('TrainingModels: ', METHOD,'-',ROBOT,'-',TASK,'-',TRAINING_SIM_REAL,'-', TRAINING_SUCCESS_FAILURE) ; ...
                  strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',char(foldname))     ;  ...
                  });
            grid on;
            grid MINOR
            plot(total_log_likelihood(sIdx,:)',COLOR(sIdx),'LineWidth',2,'Parent',subplot_1);
            axis auto
            xlim([0  length(data)]);
           %ylim([LoglikelihoodRange(1), LoglikelihoodRange(2)]);
            text(obsIdx,total_log_likelihood(sIdx,end),STATE(sIdx), 'color',COLOR(sIdx));  
            gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
            gridxy(R_State(:,2:trans),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
            hold on;
        end

        
        
        
        %{
        %Evoluting data 
        subplot_2 = subplot(2,1,2,'position', [left bottow width heigh],'Parent',gHandle_testing);
        cla(subplot_2);
        plot(data(:,1:obsIdx)','Parent',subplot_2); % for animation
        grid on;
        grid MINOR
        gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
        xlim([0  length(data)]);
        %}
        
        pause(0.000001);  
    end   
end

cd (TESTING_RESULTS_PATH);
saveas(gHandle_testing, strcat(TESTING_SIM_REAL,'_',TESTING_SUCCESS_FAILURE,'_',foldname),'jpg');
cd (rootPath);
end