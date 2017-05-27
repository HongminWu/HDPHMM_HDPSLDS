%@HongminWu March 22,2017
% testing_all_learned_models_VS_1_trial
% 4 models + 1 observation
%1. load all the learned models in a cell named as 'learnedModel'
%2. load the testing trial
function static_test_state_classification(trialID)
global  modelPath TESTING_DATASET_PATH TESTING_RESULTS_PATH PLOT_SAVE
global  METHOD ROBOT TASK STATE
global  TRAINING_SUCCESS_FAILURE TRAINING_SIM_REAL
global  TESTING_SUCCESS_FAILURE TESTING_SIM_REAL
global  THRESHOLD_PATH  TIME_STEP COLORS
global  CONST_C CONFUSION_MATRIX TIME_PERCENT

%initial parameters
trainedModel              = {};
model                     = [];
thresholdcell             = {};

expected_likelihood       = [];
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
dataLen                   = [];
model_id =  {'02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21'};
%model_id = {'23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46'};

% step-1: find out the optimal model, load the data for each state, 
%         load the expected-log_likelihood, which generated by script 'calculate_state_threshold'
STATE_STATUS        = [];
for nState = 1 : length(STATE) 
     file                      = dir([strcat(modelPath,char(STATE(nState)),'/') strcat('*',char(trialID),'.mat')]);
     optiModel                 = [];    
     optiModel                 = load(strcat(modelPath,char(STATE(nState)),'/',file(1).name));                                 %for initializating the 'meanmodel'
     trainedModel              = [trainedModel, {optiModel}];
     sData                     = data(:,R_State(nState):R_State(nState + 1));
     stateData                 = [stateData, {sData}];  
end

% step-2: calculate the cumlative log-likelihood for each state, and store in the 'stateLikelihood' array.
    stateLikelihood           = zeros(length(STATE),R_State(end));
    obsModel                  = trainedModel{1}.obsModel;
    obsModelType              = trainedModel{1}.obsModel.type;
    STATE_STATUS              =[]; 
for nState = 1:length(STATE)
    disp(strcat('Testing:',char(STATE(nState))));
    sensor                        = stateData{nState};
    data_struct.obs               = sensor;
    data_struct.true_labels       = ones(1,size(sensor,2));
    data_struct.test_cases        = 1;
    data_struct.blockSize         = ones(1,size(sensor,2));
    data_struct.blockEnd          = 1:size(sensor,2);

    if strcmp(obsModelType,'AR')
        [X,valid]                 = makeDesignMatrix(data_struct.obs,obsModel.r);
        data_struct.obs           = data_struct.obs(:,find(valid));
        data_struct.X             = X(:,find(valid));
        data_struct.blockSize     = ones(1,size(data_struct.obs,2));
        data_struct.blockEnd      = cumsum(data_struct.blockSize);
        data_struct.true_labels   = data_struct.true_labels(find(valid));
    else
        obsModel.r = 0;
    end
    state_end = [];
    state_log = [];
    for nModel = 1:length(trainedModel) %  the same data for different state model
        %testing the observed data for each learned model
        dist_struct               = trainedModel{nModel}.dist_struct; 
        theta                     = trainedModel{nModel}.theta;
        obsModel                  = trainedModel{nModel}.obsModel;
        if ~strcmp(obsModelType,'AR') 
            obsModel.r = 0;
        end
        obsModelType              = trainedModel{nModel}.obsModel.type;
        [testing_total_log_likelihood,...
            testing_neglog_c]     = observation_likelihood(data_struct,obsModelType,dist_struct,theta);   
            stateLikelihood(nModel,R_State(nState) + obsModel.r:R_State(nState + 1))     =  cumsum(testing_neglog_c); %
            state_end = [state_end, sum(testing_neglog_c)];
            state_log = [state_log; cumsum(testing_neglog_c)];
    end  
     %state classification
     [~, index] = max(state_end);
     if index == nState
         CONFUSION_MATRIX(nState, nState) = CONFUSION_MATRIX(nState, nState) + 1
     else
         CONFUSION_MATRIX(nState, index) = CONFUSION_MATRIX(nState, index) + 1
     end
    state_idx = 1:length(STATE);
    state_idx(nState) = [];
    timepos = [];
    for other_state = length(state_idx)
        diff = state_log(nState, :) - state_log( other_state,: );
        diff = abs(diff); % and then find the minimum val
        pos = find(diff==min(diff));
        timepos = [timepos, min(pos)/(R_State(nState + 1) - R_State(nState))];
    end
    timepos = min(timepos);
    TIME_PERCENT(nState) = TIME_PERCENT(nState) + timepos
end

if PLOT_SAVE
    gHandle_testing  = figure;
    %For plotting 
    
    FontSize   = 14;
    disp(strcat('Plotting! Please Wait'));
    subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
    
    fig_leg = {};
    stateLikelihood(find(stateLikelihood == 0)) = NaN;
    for sIdx = 1 : size(stateLikelihood,1)
        fig_leg{sIdx} = plot([1:length(stateLikelihood(sIdx,:))] * TIME_STEP,stateLikelihood(sIdx,:)',COLORS(sIdx),'LineWidth',3,'Parent',subplot_1);
        hold on;
    end
    title({'Cumulative Log-Likelihoods of Possible Sub-Tasks and Trained Model';...                                                                                                                        
      strcat('TrainingModels: ', METHOD,'-',ROBOT,'-',TASK,'-',TRAINING_SIM_REAL,'-', TRAINING_SUCCESS_FAILURE) ; ...
      strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',char(foldname))     ;  ...
      });
    grid on;
    axis auto
    xlabel('Time(s)','FontName','Times New Roman','FontSize',FontSize)
    ylabel('Cumulative Log-Likelihood','FontName','Times New Roman','FontSize',FontSize,'Rotation',90)
    set(gca,'xtick',[0:3:(length(data)* TIME_STEP)]);
    gridxy(R_State(:,:) * TIME_STEP,'Linestyle','--','Color',[.5 .5 .5],'linewidth',1);
    legendEntries = [fig_leg{1}, fig_leg{2}, fig_leg{3}, fig_leg{4}];
    h = legend(legendEntries,strcat('Sub-Task:',char(STATE(1))),...
                             strcat('Sub-Task:',char(STATE(2))),...
                             strcat('Sub-Task:',char(STATE(3))),...
                             strcat('Sub-Task:',char(STATE(4))),...
                             'Location','SouthWest');
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

    if (exist(TESTING_RESULTS_PATH,'dir') == 0)
        mkdir(TESTING_RESULTS_PATH);
    end
    cd (TESTING_RESULTS_PATH);
    saveas(gHandle_testing, strcat('STATE_CLASSIFICATION','_',foldname),'jpg');
end
end