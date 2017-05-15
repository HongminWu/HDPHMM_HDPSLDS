%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function static_test_state(trialID)
global rootPath modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH 
global METHOD ROBOT TASK STATE
global TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global  THRESHOLD_PATH  TIME_STEP multiModal

%initial parameters
COLOR                     = ['r', 'g', 'b','k'];
trainedModel              = {};
model                     = [];
threshold                 = [];
expected_likelihood       = [];
c                         =  5;
gHandle_testing           = figure;
sensor                    = [];
temp_log_likelihood       = zeros(length(STATE),1);
total_log_likelihood      = [];
data_struct               = struct;
left                      = 0.1; 
bottow                    = 0.1; 
width                     = 0.8; 
heigh                     = 0.4;
trans                     = 0;
stateData                 = {};
[data, R_State, foldname] = load_one_trial(trialID,TESTING_DATASET_PATH);
data                      = data_preprocessing(data); % mean = 0 and covariance = Unit matrix
stateLikelihood           = zeros(length(STATE),length(data));
for nState = 1 : length(STATE)
        %find out the optimal model for each state, this step is achieved by function "calculate_state_threshold"
        thresholdPath  = strcat(THRESHOLD_PATH, char(STATE(nState)));
        if (exist(thresholdPath,'dir') == 0)
            disp('Please define the optimal model for each state before testing? Training -> Model Selection -> Testing');
            break;
        end
        cd(thresholdPath);
        
         %load the model selection  
        model_select    = load(strcat('MODEL_SELECTION_METHOD_',char(STATE(nState)),'.mat'),'MODEL_SELECTION_METHOD');
        [~, optimalIdx] = max(model_select.MODEL_SELECTION_METHOD);
        %load the likelihoods of the optimal 
        likelihood = load(strcat(THRESHOLD_PATH, char(STATE(nState)),'/','LIKELIHOOD_',char(STATE(nState)),'.mat'),'LIKELIHOOD'); 
        state_expected_likelihood = mean(likelihood.LIKELIHOOD{optimalIdx});
        sigma                     = sqrt(var(likelihood.LIKELIHOOD{optimalIdx}));
        stateThreshold            = state_expected_likelihood - c * sigma;
        expected_likelihood       = [expected_likelihood, state_expected_likelihood];
        threshold                 = [threshold, stateThreshold];
        
        % Load Model for selected state. We will compare other observations with this.
        cd (strcat(modelPath,char(STATE(nState))));
        file           = dir('*.mat');
        optiModel      = [];    
        optiModel      = load(file(optimalIdx).name);                                 %for initializating the 'meanmodel'
        trainedModel   = [trainedModel, {optiModel}];
        sData          = data(:,R_State(nState):R_State(nState + 1));
        stateData      = [stateData, {sData}];
end

        obsModel       = trainedModel{1}.obsModel;
        obsModelType   = trainedModel{1}.obsModel.type;
for nState = 1:length(STATE)
    disp(strcat('Testing:',char(STATE(nState))));
    sensor                  = stateData{nState};
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
    
    for nModel = 1:length(trainedModel)
        %testing the observed data for each learned model
        dist_struct         = trainedModel{nModel}.dist_struct; 
        theta               = trainedModel{nModel}.theta;
        obsModel            = trainedModel{nModel}.obsModel;
        obsModelType        = trainedModel{nModel}.obsModel.type;
        [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);       
        stateLikelihood(nModel,R_State(nState) + obsModel.r:R_State(nState + 1))     =  cumsum(testing_neglog_c); %
    end
end

   %For plotting  
    FontSize   = 9;
    disp(strcat('Plotting! Please Wait'));
    subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
    fig_expected_likelihood = plot([1:length(expected_likelihood)] * TIME_STEP, expected_likelihood','m--','LineWidth',2,'Parent',subplot_1);
    hold on;
    fig_threshold = plot([1:length(threshold)] * TIME_STEP, threshold','c--','LineWidth',2,'Parent',subplot_1);     
    hold on;

    fig_leg = {};
    for sIdx = 1 : size(stateLikelihood,1)
        title({'Cumulative Log-Likelihoods of Possible Sub-Tasks and Trained Model';...                                                                                                                        
              strcat('TrainingModels: ', METHOD,'-',ROBOT,'-',TASK,'-',TRAINING_SIM_REAL,'-', TRAINING_SUCCESS_FAILURE) ; ...
              strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',char(foldname))     ;  ...
              });
        grid on;
        grid MINOR
        fig_leg{sIdx} = plot([1:length(stateLikelihood(sIdx,:))] * TIME_STEP,stateLikelihood(sIdx,:)',COLOR(sIdx),'LineWidth',2,'Parent',subplot_1);
        axis auto
        xlabel('Time(s)','FontName','Times New Roman','FontSize',FontSize)
        ylabel('Cumulative Log-Likelihood','FontName','Times New Roman','FontSize',FontSize,'Rotation',90)
        set(gca,'xtick',[0:3:(length(data)* TIME_STEP)]);
        gridxy(R_State(:,:) * TIME_STEP,'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
        hold on;
    end

    legendEntries = [fig_expected_likelihood, fig_threshold, fig_leg{1}, fig_leg{2}, fig_leg{3}, fig_leg{4}];
    h = legend(legendEntries,'Expected-Log\_likelihood of Each Sub-Task','Threshold of Each Sub-Task','Sub-Task:Approach','Sub-Task:Rotation','Sub-Task:Insertion','Sub-Task:Mating', 'Location','SouthEast');
    set(h,'Fontsize',FontSize);
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
%set(gHandle_testing,'PaperUnits','centimeters','PaperPosition',[14, 19, 10, 15],'PaperPositionMode','manual');
if (exist(TESTING_RESULTS_PATH,'dir') == 0)
    mkdir(TESTING_RESULTS_PATH);
end
set(gHandle_testing, 'position', [0,0,1800,1000]);
cd (TESTING_RESULTS_PATH);
saveas(gHandle_testing, strcat('_',foldname),'jpg');
cd (rootPath);
end