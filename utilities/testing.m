%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function testing(trialID)
global rootPath modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH 
global STATE TRUNCATION_STATES TRUNCATION_COMPONENTS
global TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global COLOR TRAINING_PLOT_ON TESTING_PLOT_ON
%% Training Model
learnedModel = {};
model = [];
for nlearnedModel = 1 : length(STATE)
    % Load Model for selected state. We will compare other observations with this.
        cd (strcat(modelPath,char(STATE(nlearnedModel))));
        file         = dir('*.mat');
        meanmodel = [];
        meanmodel    = load(file(2).name);                                 %for initializating the 'meanmodel'
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
 %      mean_loglike = plot_traning_data(training_model_path, TRAINING_PLOT_ON);
end

[data, R_State, foldname] = load_one_trial(trialID,TESTING_DATASET_PATH);
cd(rootPath);
%% testing
gHandle_testing = figure;
title(strcat('Testing: The likelihood of testing trial w.r.t the learned model'));
sensor = [];
temp_log_likelihood = zeros(length(STATE),1);
total_log_likelihood = [];
%% Inference using training model and testing trial
settings.Kz             =  TRUNCATION_STATES;  % nStats
settings.Ks             =  TRUNCATION_COMPONENTS;  % nComponents of mixture gaussian
data_struct             = struct;
left=0.1; bottow=0.1; width=0.8; heigh=0.4;
trans = 1;
plot_data = data_preprocessing(data);
for obsIdx = 1:length(data)   
    disp(strcat('Testing: obsIdx = ',num2str(obsIdx),'/',num2str(length(data))));
    %load different learned model at different state
    if ismember(obsIdx, R_State(:,2:end-1))
        sensor = [];
        gridxy(R_State(:,2:end-1),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
        trans = trans + 1;
    end
    sensor                  = [sensor, data(:,obsIdx)];
    [sensor, meanSigma]     = data_preprocessing(sensor); % mean = 0 and covariance = Unit matrix
  
    data_struct.obs         = sensor;
    data_struct.true_labels = ones(1,size(sensor,2));
    data_struct.test_cases  = 1;
    data_struct.blockSize   = ones(1,size(sensor,2));
    data_struct.blockEnd    = 1:size(sensor,2);
    %testing the observed data for each learned model
    for nlearnedModel       = 1 : length(STATE)
        dist_struct         = learnedModel{nlearnedModel}.dist_struct; 
        theta               = learnedModel{nlearnedModel}.theta;
        obsModel            = learnedModel{nlearnedModel}.obsModel;
        obsModelType        = learnedModel{nlearnedModel}.obsModel.type;
        dist_struct.pi_init = [1.0,zeros(1,TRUNCATION_STATES-1)];
        [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);

%          [likelihood log_normalizer] = compute_likelihood(data_struct,theta,obsModelType,settings.Kz,settings.Ks);
%          testing_total_log_likelihood = sum(squeeze(log_normalizer));
%         testing_total_log_likelihood = exp(max(squeeze (max(max(log_normalizer,[],1),[],2)))); 
        
        temp_log_likelihood(nlearnedModel,1) = testing_total_log_likelihood ; %
    end
    total_log_likelihood = [total_log_likelihood, temp_log_likelihood]; %store the results
    
    if TESTING_PLOT_ON            
        %info: subplot('position',[left bottom width height]) range from 0.0 to 1.0
        clf(gHandle_testing);
        subplot_1 = subplot(2,1,1,'position', [left bottow + heigh + 0.02 width heigh],'Parent',gHandle_testing);
        for sIdx = 1 : size(total_log_likelihood,1)
            title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...
                  strcat('TrainingModels: ', TRAINING_SIM_REAL,'\_', TRAINING_SUCCESS_FAILURE); ...
                  strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',foldname)});
            grid on;
            grid MINOR
            plot(total_log_likelihood(sIdx,:)',COLOR(sIdx),'LineWidth',2,'Parent',subplot_1);
            xlim([0  length(data)]);
           % ylim([LoglikelihoodRange(1), LoglikelihoodRange(2)]);
            text(obsIdx,total_log_likelihood(sIdx,end),STATE(sIdx), 'color',COLOR(sIdx));  
            gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
            gridxy(R_State(:,2:trans),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
            hold on;
        end
        
        subplot_2 = subplot(2,1,2,'position', [left bottow width heigh],'Parent',gHandle_testing);
        plot(plot_data(:,1:obsIdx)','Parent',subplot_2); % for animation
        %plot(data','Parent',subplot_2);
        grid on;
        grid MINOR
        gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
        xlim([0  length(data)]);
        pause(0.000001);
    elseif ismember(obsIdx, R_State(:,end)) %the ending point
        subplot(2,1,1,'position', [left bottow + heigh  + 0.02 width heigh]);
         for sIdx = 1 : size(total_log_likelihood,1)
            title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...
                  strcat('TrainingModels: ', TRAINING_SIM_REAL,'\_', TRAINING_SUCCESS_FAILURE,'\_',foldname); ...
                  strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',foldname)});
            grid on;
            grid MINOR
            xlim([0  length(data)]);
            plot(total_log_likelihood(sIdx,:)',COLOR(sIdx),'LineWidth',1);
            text(size(total_log_likelihood,2),total_log_likelihood(sIdx,end),STATE(sIdx), 'color',COLOR(sIdx));  
            gridxy(R_State(:,2:trans),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
            hold on;
         end
        subplot(2,1,2,'position', [left bottow width heigh]);
        plot(data');
        grid on;
        grid MINOR
        xlim([0  length(data)]);
    end   
end

cd (TESTING_RESULTS_PATH);
saveas(gHandle_testing, strcat(TESTING_SIM_REAL,'_',TESTING_SUCCESS_FAILURE,'_',foldname),'jpg');
cd (rootPath);
end
function mean_loglike = plot_traning_data(training_model_path, trainPlotON)
cd (training_model_path);
file = dir('*.mat');

input_data = {};
align_data = [];
nbData = [];
temp_data = load(file(1).name);
temp_loglike = squeeze(cumsum(temp_data.output_neglog_c{:}));
L = length(temp_loglike);
for i=1:length(file) % a folder per demonsration
    training_data = load(file(i).name);
    %training_loglike = squeeze(cumsum(training_data.output_neglog_c{:}));
    training_loglike = squeeze(training_data.output_neglog_c{:});
    training_loglike(L+1:end,:) = [];
    input_data{i} = training_loglike;
    nbData = [nbData, length(input_data{i})]; % so as to find out the max length
end
for n=1:length(input_data)
    align = spline(1:size(input_data{n},2), input_data{n}, linspace(1,size(input_data{n},2), max(nbData)));
    align_data = [align_data; align];%Resampling alignment
    if trainPlotON 
        plot(align');
    end
end

mean_loglike = mean(align_data);
if length(input_data)==1
    mean_loglike = align_data;
end
if trainPlotON
    plot(mean_loglike','b','LineWidth',3);
end
end
