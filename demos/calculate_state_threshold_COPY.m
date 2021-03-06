%step1: Choose N trials for training, after this step, we can obtain N trained models
%step2: Choose 1 testing trial w.r.t N-1 trained model(except 1 testing trial).  
        %After this step, we can get (N - 1) log_likelihood vectors for each trained model
%step3: Calculate the mean vector and variance vector of log_likelihood for each trained model. 
        %After this step, we can get the 1 mean vector and 1 variance vector for each trained model
%step4: Choose the trained model by compare the mean value of ending point of each mean vector: max(meanVector(end))
%step5: the expected log_likelihood is generated by trained model that has the maximum mean value of the ending point through the same file_id trial
        %as the testing trial
%step6: threshold = mean - c*standard variation
function calculate_state_threshold_COPY()
clc;
close all;
clear global;
global_variables;
global modelPath STATE TRAINING_DATASET_PATH THRESHOLD_PATH DIS_PERIOD
[APPROACH, ROTATION, INSERTION, MATING] = deal(1,2,3,4);
c       = 0.5; %for calculating the threshood
for nState = APPROACH : APPROACH
    data_struct = struct;
    % Load Model for selected state.
    learnedModel = {};
    FILE_ID      = [];
    cd (strcat(modelPath,char(STATE(nState))));
    file            = dir('*.mat');
    %load all the trained models and find the trial from the learned model
    for nlearnedModel = 1 : length(file)
        fileID       = uint16(str2double(file(nlearnedModel).name((length(file(nlearnedModel).name)-5):(length(file(nlearnedModel).name)-4)))); %find the file id
        FILE_ID      = [FILE_ID, fileID];
        model    = load(file(nlearnedModel).name); 
        learnedModel = [learnedModel, {model}];
    end
    
    %cd saving path
    thresholdPath  = strcat(THRESHOLD_PATH, char(STATE(nState)));
    if (exist(thresholdPath,'dir') == 0)
        mkdir(thresholdPath);
    end
    cd(thresholdPath);
    delete *;
    
    %left one cross validataion
    MODEL_SELECTION_METHOD = [];
    LIKELIHOOD     = {};
    for n = 1:length(FILE_ID) % for each testing file 
        modelIdx = 1:length(learnedModel);
        modelIdx(n) = [];
        
        %step2-1: load 1 testing trail
        [sensor, R_State, foldname] = load_one_trial(FILE_ID(n),TRAINING_DATASET_PATH);
        [data, meanSigma] = data_preprocessing(sensor);
        data = data(:,R_State(nState):R_State(nState + 1));
        
        %step2-2: load N-1 learned models
        log_likelihood = zeros(length(modelIdx),length(data));
        for m = 1:length(modelIdx) % for each trained model
           % calculate the likelihoood
           obsModel            = learnedModel{modelIdx(m)}.obsModel;
           obsModelType        = learnedModel{modelIdx(m)}.obsModel.type;
           r = 1;
           if strcmp(obsModelType,'AR')
               r = r + obsModel.r;
           end
           sensor = [];
           for obsIdx = r:length(data)
               if rem(obsIdx, DIS_PERIOD) == 0
                  disp(strcat(char(STATE(nState)),': ',num2str(n),'/',num2str(length(FILE_ID)),'-Model-',num2str(m),'/',num2str(length(modelIdx)),'-obsIdx-',num2str(obsIdx),'/',num2str(length(data))));
               end
               if strcmp(obsModelType,'AR') & obsIdx == r
                    sensor =  data(:, (obsIdx - obsModel.r): obsIdx);
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
                dist_struct         = learnedModel{modelIdx(m)}.dist_struct; 
                theta               = learnedModel{modelIdx(m)}.theta;
                obsModel            = learnedModel{modelIdx(m)}.obsModel;
                obsModelType        = learnedModel{modelIdx(m)}.obsModel.type;
                [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);       
                log_likelihood(m,obsIdx) = testing_total_log_likelihood ; %
           end  
        end   %.w.r.t N-1 models  
        LIKELIHOOD = [LIKELIHOOD, {log_likelihood}];
        
        %calculate the mean likelihood
       %#-1
       mean_likelihood        = mean(log_likelihood);
       MODEL_SELECTION_METHOD = [MODEL_SELECTION_METHOD, mean_likelihood(end)];
       %#-2
        %MODEL_SELECTION_METHOD = [MODEL_SELECTION_METHOD, max(log_likelihood(:,end))];
    end  % end for one state 
    
    cd(thresholdPath);
    %save the MEAN_LIKELIHOOD_END for selecting the optimal model
    fileName = strcat('MODEL_SELECTION_METHOD_',char(STATE(nState)));
    save(fileName,'MODEL_SELECTION_METHOD');
    
    %save all the likelihood for each trained model
    fileName = strcat('LIKELIHOOD_',char(STATE(nState)));
    save(fileName,'LIKELIHOOD');
    
    [~,optimalIndex] = max(MODEL_SELECTION_METHOD);  %%%%%%%find the optimal model 
    
    %calculate the the threshold
    optiLIKELIHOOD = LIKELIHOOD{optimalIndex};
    mu             = mean(optiLIKELIHOOD);     
    sigma          = sqrt(var(optiLIKELIHOOD));
    threshold      = mu - c * sigma;
    %save threshold
    fileName = strcat('THRESHOLD_',char(STATE(nState)));
    save(fileName,'threshold');
    
    expected_likelihood = mu; %Expectation Maximization
    
    cd(thresholdPath);
            
    %save the expected likelihood
    fileName = strcat('EXPECTED_LIKELIHOOD_',char(STATE(nState)));
    save(fileName,'expected_likelihood');
    
    gh = figure;
    plot(threshold','r--');
    hold on;
    plot(expected_likelihood','b','LineWidth',2);
    title(char(STATE(nState)));
    legend('Threshold','Expected log-likelihood','Location','NorthWest')
    saveas(gh, char(STATE(nState)),'jpg');
end
end