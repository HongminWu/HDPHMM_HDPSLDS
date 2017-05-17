%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function test_state_COPY(trialID)
global rootPath modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH 
global METHOD ROBOT TASK STATE TRUNCATION_STATES TRUNCATION_COMPONENTS
global TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global TESTING_PLOT_ON THRESHOLD_PATH DIS_PERIOD TIME_STEP
                
COLOR                  = ['r', 'g', 'b','k'];
%% Training Model
trainedModel           = {};
model                  = [];
threshold              = [];
expected_likelihood    = [];
c                      =  5; 
for nState = 1 : length(STATE)
    
        %find out the optimal model for each state, this step is achieved by function "calculate_state_threshold"
        thresholdPath  = strcat(THRESHOLD_PATH, char(STATE(nState)));
        if (exist(thresholdPath,'dir') == 0)
            disp('Please define the optimal model for each state before testing? Training -> Model Selection -> Testing');
            break;
        end
        
        cd(thresholdPath);
        % load the model selection
        model_select = load(strcat('MODEL_SELECTION_METHOD_', chat(STATE(nState)), '.mat'));
        [~, optimalIdx] = max(model_select.MODEL_SELECTION_METHOD);
        %load the optimal likelihoods  
        likelihood                = load(strcat(THRESHOLD_PATH, char(STATE(nState)),'/','LIKELIHOOD_',char(STATE(nState)),'.mat')); 
        state_expected_likelihood = mean(likelihood.LIKELIHOOD{optimalIdx});
        sigma                     = sqrt(var(likelihood.LIKELIHOOD{optimalIdx}));
        stateThreshold            = state_expected_likelihood - c * sigma;
        
        % Load Model for selected state. We will compare other observations with this.
        cd (strcat(modelPath,char(STATE(nState))));
        file         = dir('*.mat');
        optiModel = [];    
        optiModel    = load(file(optimalIdx).name);                                 %for initializating the 'meanmodel'
        trainedModel                   = [trainedModel, {optiModel}];
end
[data, R_State, foldname] = load_one_trial(trialID,TESTING_DATASET_PATH);
cd(rootPath);
%% testing
gHandle_testing = figure;
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

obsModel            = trainedModel{nState}.obsModel;
obsModelType        = trainedModel{nState}.obsModel.type;
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
        trans     = trans + 1;
        if strcmp(obsModelType,'AR')
            sensor =  data(:, (obsIdx - obsModel.r): obsIdx);
        end  
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
    for nState       = 1 : length(STATE)
        dist_struct         = trainedModel{nState}.dist_struct; 
        theta               = trainedModel{nState}.theta;
        obsModel            = trainedModel{nState}.obsModel;
        obsModelType        = trainedModel{nState}.obsModel.type;
        [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);       
        temp_log_likelihood(nState,1) = testing_total_log_likelihood ; %
    end
    total_log_likelihood = [total_log_likelihood, temp_log_likelihood]; %store the results
    
    if TESTING_PLOT_ON & ((rem(obsIdx, DIS_PERIOD) == 0) | ismember(obsIdx, R_State(:,end)))  
        
        clf(gHandle_testing);
        %subplot_1 = subplot(2,1,1,'position', [left bottow + heigh  width heigh],'Parent',gHandle_testing);
        subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
        fig_expected_likelihood = plot([1:length(expected_likelihood)] * TIME_STEP, expected_likelihood','m--','LineWidth',2,'Parent',subplot_1);
        hold on;
        fig_threshold = plot([1:length(threshold)] * TIME_STEP, threshold','c--','LineWidth',2,'Parent',subplot_1);     
        hold on;
        
        fig_leg = {};
        for sIdx = 1 : size(total_log_likelihood,1)
            title({'Cumulative Log-Likelihoods of Possible Sub-Tasks and Trained Model';...                                                                                                                        
                  strcat('TrainingModels: ', METHOD,'-',ROBOT,'-',TASK,'-',TRAINING_SIM_REAL,'-', TRAINING_SUCCESS_FAILURE) ; ...
                  strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',char(foldname))     ;  ...
                  });
            grid on;
            grid MINOR
            fig_leg{sIdx} = plot([1:obsIdx - 1] * TIME_STEP,total_log_likelihood(sIdx,:)',COLOR(sIdx),'LineWidth',2,'Parent',subplot_1);
            axis auto
%           xlim([0  length(data)]* TIME_STEP);
            xlabel('Time(s)','FontName','Times New Roman','FontSize',14)
            ylabel('Cumulative Log-Likelihood','FontName','Times New Roman','FontSize',14,'Rotation',90)
            set(gca,'xtick',[0:3:(length(data)* TIME_STEP)]);
            
%           text(obsIdx,total_log_likelihood(sIdx,end),STATE(sIdx), 'color',COLOR(sIdx));  
            gridxy(obsIdx * TIME_STEP,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
            gridxy(R_State(:,1:trans) * TIME_STEP,'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
            hold on;
        end

        legendEntries = [fig_expected_likelihood, fig_threshold, fig_leg{1}, fig_leg{2}, fig_leg{3}, fig_leg{4}];
        h = legend(legendEntries,'Expected-Log\_likelihood of Each Sub-Task','Threshold of Each Sub-Task','Sub-Task:Approach','Sub-Task:Rotation','Sub-Task:Insertion','Sub-Task:Mating', 'Location','SouthEast');
        set(h,'Fontsize',16);
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
set(gHandle_testing, 'unit', 'normalized', 'position', [0,0,1,1]);
if (exist(TESTING_RESULTS_PATH,'dir') == 0)
    mkdir(TESTING_RESULTS_PATH);
end
cd (TESTING_RESULTS_PATH);
saveas(gHandle_testing, strcat('_',foldname),'jpg');
cd (rootPath);
end