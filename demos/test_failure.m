function test_failure()
    global_variables;
    cd (strcat(modelPath,'APPROACH'));
    trialID                     = 1;
    TESTING_DATASET             = 'REAL_HIRO_ONE_SA_ERROR_CHARAC'; %SIM_HIRO_ONE_SA_FAILURE, REAL_HIRO_ONE_SA_ERROR_CHARAC
    TESTING_DATASET_PATH        = strcat(rootDATApath,TESTING_DATASET);
    file                        = dir('*.mat'); 
    trained_model_name          = file(1).name;
    learnedModel                = load(trained_model_name);  
    [data, R_State, foldname]   = load_one_trial(trialID,TESTING_DATASET_PATH);
    title(strcat('Testing: The likelihood of testing trial w.r.t the learned model'));
    gHandle_testing             = figure;
    sensor                      = [];
    total_log_likelihood        = [];
    data_struct                 = struct;
    data                        = data_preprocessing(data); % mean = 0 and covariance = Unit matrix
    for obsIdx = 1 : length(data) 
        disp(strcat('Testing: obsIdx = ',num2str(obsIdx),'/',num2str(length(data))));
        sensor                  = [sensor, data(:,obsIdx)];
        data_struct.obs         = sensor;
        data_struct.true_labels = ones(1,size(sensor,2));
        data_struct.test_cases  = 1;
        data_struct.blockSize   = ones(1,size(sensor,2));
        data_struct.blockEnd    = 1:size(sensor,2);
        dist_struct             = learnedModel.dist_struct; 
        theta                   = learnedModel.theta;
        obsModel                = learnedModel.obsModel;
        obsModelType            = learnedModel.obsModel.type;
        [testing_total_log_likelihood,
            testing_neglog_c]   = observation_likelihood(data_struct,obsModelType,dist_struct,theta);       
        total_log_likelihood    = [total_log_likelihood; testing_total_log_likelihood]; %

        %plot
        clf(gHandle_testing);
        subplot_1 = subplot(1,1,1,'Parent',gHandle_testing);
        title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...
              strcat('TrainingModels: ', TRAINING_SIM_REAL,'\_', TRAINING_SUCCESS_FAILURE,'\_',trained_model_name); ...
              strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',foldname)});
        grid on;
        grid MINOR
        plot(total_log_likelihood,'b-','LineWidth',2,'Parent',subplot_1);
        axis auto
        xlim([0  length(data)]);
        gridxy(obsIdx,'Linestyle','-','Color',[.5 .5 .5],'linewidth',0.6);
        hold on;
        pause(0.000001);    
%         if ismember(obsIdx, R_State(:,end)) %the ending point
%             subplot(2,1,1,'position', [left bottow width heigh]);
%             title({'Given the learned model for each state, calculate the log-likelihood of the continuous observations';...
%                   strcat('TrainingModels: ', TRAINING_SIM_REAL,'\_', TRAINING_SUCCESS_FAILURE,'\_',trained_model_name); ...
%                   strcat('TestingTask: '   , TESTING_SIM_REAL ,'\_', TESTING_SUCCESS_FAILURE ,'\_',foldname)});
%             grid on;
%             grid MINOR
%             xlim([0  length(data)]);
%             plot(total_log_likelihood(sIdx,:)',COLOR(sIdx),'LineWidth',1);
%             text(size(total_log_likelihood,2),total_log_likelihood(sIdx,end),STATE(sIdx), 'color',COLOR(sIdx));  
%             gridxy(R_State(:,2:trans),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);
%             hold on;
%         end
    end
end