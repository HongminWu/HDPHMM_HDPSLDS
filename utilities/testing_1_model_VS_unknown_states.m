%@HongminWu 02 23,2017
% This script is used for specific state recognition by comparing the likelihood 
%1. caculating the average likelihood of all training trails
%2. 
function testing_1_model_VS_unknown_states()
close all;
global STATE
%% declaration
stateName = {'approach', 'rotation', 'insection', 'mating'};
[approach, rotation, insertion, mating] = deal(1, 2, 3, 4);
color = ['r', 'g', 'b','c'];
figFORMAT = 'jpg';
% Select training model you want to use to compare with other observations.
training_state = rotation; 
% Paths
rootPath             = 'F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox\HDPHMM_HDPSLDS';
training_model_path  = strcat(rootPath,'\output_sim_success\success_',char(stateName(training_state)));
testing_dataset_path = strcat(rootPath,'\logs\SIM_HIRO_ONE_SA_SUCCESS');
testing_results_path = strcat(rootPath,'\testing_results');

%% Training Model
figure;
title(strcat('Training: The average likelihood of\_ ',char(stateName(training_state))));
hold on;
% Load Model for selected state. We will compare other observations with this.
load(strcat(training_model_path,'\success_output_10000_1.mat'));    % load the parameters of learned HMM model
mean_loglike = plot_traning_data(training_model_path);

%% loading data
%Get Observation data from testing trials
id = 10 + 1; %testing trail
cd(testing_dataset_path); %test dataset
folders = dir;     
dataType = {'R_Torques'};
for i=1:size(id, 2) % a folder per demonsration
    if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git') && ~strcmp(folders(id(i)).name,'ERROR'))
        foldname = folders(id(i)).name;
        cd(foldname);
        data = [];
        file = dir([char(dataType(1)) '.dat*']);
        filename = strcat('/',foldname,'/',file(1).name);
        raw_data = load(filename);
        d = raw_data(:,2:end);  %delete the time column   
        if strcmp(dataType(1), 'R_Torques')             
            d = unique(d,'rows','stable');
            d = [d,[d(1,:);diff(d)]];
        end

        data = [data, d'];
        R_State = [];
        file = dir([char('R_State') '.dat*']);
        filename = strcat('/',foldname,'/',file(1).name);
        tState = load(filename);
        for in = 1:length(tState)
            idex= find(raw_data(:,1) == tState(in));
            R_State = [R_State,idex];
        end
        R_State = [R_State,length(data)];
    end
end

%% testing
gHandle_testing = figure;
title(strcat('Testing: The likelihoods of unknown states w.r.t the learned model (',char(stateName(training_state)),')'));
for testing_state = 1:length(stateName)   % compare states
    disp(strcat('Training:',char(stateName(training_state)),'--Testing:',char(stateName(testing_state)),'(color:',color(testing_state),')'));
    %data preprocessing
    %split the specific state    
    sensor = data(:,R_State(testing_state):R_State(testing_state+1));
    % Smooth trajextory
    radius = 20;
    [dims,len] = size(sensor);
    smoothed = zeros(dims,len);
    for j=1:dims
         for k=1:len
             low = max(1,k-radius);
             high = min(len,k+radius);
             smoothed(j,k) = mean(sensor(j,low:high));
         end
    end
     sensor = smoothed;
    % Adjust each dim to mean 0
    mY = mean((sensor'));
    for j=1:length(mY)
        sensor(j,:) = sensor(j,:) - mY(j);
    end
    %Renormalize so for each feature, the variance of the first diff is 1.0
    vY = var(diff(sensor'));
    for j=1:length(vY)
        sensor(j,:) = sensor(j,:) ./ sqrt(vY(j));
    end
    meanSigma = 5.0 * cov(diff(sensor'));  %If bad segmentation, try values between 0.75 and 5.0
    for i=1:size(meanSigma,1)
        for j=1:size(meanSigma,2)
            if(i~=j)
                meanSigma(i,j) = 0;
            end
        end
    end
    
%% Inference using training model and testing trial
    settings.Kz             = 10;  % nStats
    settings.Ks             = 10;  % nComponents of mixture gaussian
    
    data_struct             = struct;
    data_struct.obs         = sensor;
%     data_struct.X           = output_data_struct{:}.X(:,1:length(sensor)-1) ;
    data_struct.true_labels = ones(1,size(sensor,2));
    data_struct.test_cases  = 1;
    data_struct.blockSize   = ones(1,size(sensor,2));
    data_struct.blockEnd    = 1:size(sensor,2); 
    dist_struct             = output_dist_struct{:}; 
    theta                   = output_theta{:};
    obsModel                = output_obsModel{:};
    obsModelType            = output_obsModel{:}.type;
    
    
   [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);
    hold on;
    plot(squeeze(cumsum(testing_neglog_c)),color(testing_state),'LineWidth',3);
%     plot(squeeze(testing_neglog_c),color(testing_state),'LineWidth',3);
    pause(2);
end
cd (testing_results_path);
picName = strcat(char(stateName(testing_state)),'_',char(stateName(training_state)));
saveas(gHandle_testing, picName,figFORMAT);
cd (rootPath);
end
function mean_loglike = plot_traning_data(training_model_path)
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
    training_loglike = squeeze(cumsum(training_data.output_neglog_c{:}));
%     training_loglike = squeeze(training_data.output_neglog_c{:});
    training_loglike(L+1:end,:) = [];
    input_data{i} = training_loglike;
    nbData = [nbData, length(input_data{i})]; % so as to find out the max length
end
for n=1:length(input_data)
    align = spline(1:size(input_data{n},2), input_data{n}, linspace(1,size(input_data{n},2), max(nbData)));
    align_data = [align_data; align];%Resampling alignment
    plot(align');
end

mean_loglike = mean(align_data);
if length(input_data)==1
    mean_loglike = align_data;
end
plot(mean_loglike','b','LineWidth',3);
end