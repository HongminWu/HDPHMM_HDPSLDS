%@HongminWu 02 23,2017
function demo_testing()
close all;
clear;
clc;

%% Parameters
[approach, rotation, insertion, mating] = deal(1, 2, 3, 4);
color = ['r', 'g', 'k','c'];

% Paths
training_model_path='F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox\HDPHMM_HDPSLDS\output_sim_success\success_';

%% Training Model
% Select training model you want to use to compare with other observations.
training_state = insertion; 

switch training_state
    case approach
        training_seg_name = 'approach';
    case rotation
        training_seg_name = 'rotation';
    case insertion
        training_seg_name = 'insection';
    case mating
        training_seg_name = 'mating';
end

% Load Model for selected state. We will compare other observations with this.
tp=strcat(training_model_path,training_seg_name,'\success_output_1000_1.mat');
load(tp);    % load the parameters of learned HMM model
mean_loglike = plot_traning_data(strcat(training_model_path,training_seg_name));

%% loading data
compare_loglike_normalizer = {};
for testing_state = approach:mating   % compare states
    iter = 50;
    switch testing_state
        case approach
            testing_seg_name = 'approach';
        case rotation
            testing_seg_name = 'rotation';
        case insertion
            testing_seg_name = 'insection';
        case mating
            testing_seg_name = 'mating';
    end
    
    
    %% Get Observation data from testing trials
    cd('F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox\HDPHMM_HDPSLDS\logs\SIM_HIRO_ONE_SA_SUCCESS'); %test dataset
    DataCell = {};
    R_State = {};
    folders = dir;
    
    id = 3 + 6; %testing trail
    for i=1:size(id, 2) % a folder per demonsration
        if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git') && ~strcmp(folders(id(i)).name,'ERROR'))
            foldname = folders(id(i)).name;
            cd(foldname);
            data = [];
            file = dir(['R_Torques' '.dat']);
            filename = strcat('/',foldname,'/',file(1).name);
            raw_data = load(filename);
            d = raw_data(:,2:end);  %delete the time column
%             d = smoothWrenchFilter(d); % wrenchVec = d;
            d = unique(d,'rows','stable');
            d = [d,[d(1,:);diff(d)]]; % first-order derivative
            data = [data, d];
            DataCell = [DataCell; {data}];
            
            iState = [];
            file = dir([char('R_State') '.dat*']);
            filename = strcat('/',foldname,'/',file(1).name);
            tState = load(filename);
            for in = 1:length(tState)
                idex= find(raw_data(:,1) == tState(in));
                iState = [iState,idex];
            end
            iState = [iState,size(data,1)];
            R_State = [R_State; {iState}];
            cd ..
        end
    end
    
    % Get length of observations for testing based on training state
%    for i = 1: length(DataCell) % split the specific state    
     training_length=length( R_State{1}(training_state):R_State{1}(training_state + 1) );
     DataCell{1} = DataCell{1}(R_State{1}(testing_state):R_State{1}(testing_state + 1),:);
     if(length(DataCell)>training_length)
         DataCell=DataCell(1:training_length);
     elseif(length(DataCell)<training_length)
         training_length=length(DataCell{1});
     end
 %   end
    
    % alignment
    sensor = [];
    nbData = [];
    for i = 1: length(DataCell)
        DataCell{i} = DataCell{i}';
        nbData = [nbData, length(DataCell{i})]; % so as to find out the max length
    end
    for n=1:length(DataCell)
        %DataCell{n} = spline(1:size(DataCell{n},2), DataCell{n}, linspace(1,size(DataCell{n},2), max(nbData))); %Resampling alignment
        %         DataCell{n} = spline(1:size(DataCell{n},2), DataCell{n}, linspace(1,size(DataCell{n},2), length(mean_loglike))); %Resampling alignment
        sensor = [sensor; DataCell{n}];
    end
    
    %Smooth trajextory
    %         radius = 20;
    %         [dims,len] = size(sensor);
    %         smoothed = zeros(dims,len);
    %         for j=1:dims
    %             for k=1:len
    %                 low = max(1,k-radius);
    %                 high = min(len,k+radius);
    %                 smoothed(j,k) = mean(sensor(j,low:high));
    %             end
    %         end
    %         sensor = smoothed;
    
    %% Adjust each dim to mean 0
    mY = mean((sensor'));
    for j=1:length(mY)
        sensor(j,:) = sensor(j,:) - mY(j);
    end
    %Renormalize so for each feature, the variance of the first diff is 1.0
    vY = var(diff(sensor'));
    for j=1:length(vY)
        sensor(j,:) = sensor(j,:) ./ sqrt(vY(j));
    end
    meanSigma = 1.0 * cov(diff(sensor'));  %If bad segmentation, try values between 0.75 and 5.0
    sensor = diff(sensor')'; %let the data mean=0, variance = 1
    for i=1:size(meanSigma,1)
        for j=1:size(meanSigma,2)
            if(i~=j)
                meanSigma(i,j) = 0;
            end
        end
    end
    
    %% Inference using training model and testing trial
    settings.Kz             = 10;  % nStats
    settings.Ks             =  1;  % nComponents of mixture gaussian
    
    data_struct             = struct;
    data_struct.obs         = sensor;
    data_struct.true_labels = ones(1,length(sensor));
    data_struct.test_cases  = 1;
    data_struct.blockSize   = ones(1,length(sensor));
    data_struct.blockEnd    = 1:length(sensor);
    
    dist_struct             = output_dist_struct{:};
    
    theta                   = output_theta{:};
    obsModel                = output_obsModel{:};
    obsModelType            = output_obsModel{:}.type;
    
    %gHandle_obs = figure;
    %while(iter > 0)
        disp(strcat('Training State: ', training_state, 'Iter: ', num2str(501 - iter),'/',testing_seg_name));
        
        % Block sample (z_{1:T},s_{1:T})|y_{1:T}
        [stateSeq INDS stateCounts] = sample_zs(data_struct,dist_struct,theta,obsModelType);
        
        % Create sufficient statistics:
        Ustats = update_Ustats(data_struct,INDS,stateCounts,obsModelType);
        
        % Sample theta
        % theta = sample_theta(theta,Ustats,obsModel);
        
        [testing_total_log_likelihood, testing_neglog_c] = observation_likelihood(data_struct,obsModelType,dist_struct,theta);
        % [testing_likelihood testing_loglike_normalizer] = compute_likelihood(data_struct,theta,obsModelType,settings.Kz,settings.Ks);
        
        %         figure(gHandle_obs);
        %         plot(mean_loglike,'b','LineWidth',2);
        %         hold on;
        %         plot(squeeze(cumsum(testing_neglog_c)),'b-','LineWidth',1);
        %         title('Testing:Cumulative Likelihood');
        % drawnow;
        iter = iter -1;
        
    %end
    plot(squeeze(cumsum(testing_neglog_c)),color(testing_state),'LineWidth',3);
    hold on;
    pause
end
% plot(mean_loglike,'b','LineWidth',2);
hold on;

%     legend('training\_mean\_loglike\_of\_Rotation','testing\_loglike\_of\_approach','testing\_loglike\_of\_rotation'...
%            ,'testing\_loglike\_of\_insection','testing\_loglike\_of\_mating','Location','SouthWest');

cd ('../../testing_results');
%     saveas(gHandle_obs, foldname,'jpg');
cd ..
end
function mean_loglike = plot_traning_data(success_path)
% cd (success_path);
file = dir(success_path);
file(1:3)=[];

input_data = {};
align_data = [];
nbData = [];
temp_data = load(strcat(success_path,'/',file(1).name));
temp_loglike = squeeze(cumsum(temp_data.output_neglog_c{:}));
L = length(temp_loglike);
for i=1:length(file) % a folder per demonsration
    training_data = load(file(i).name);
    training_loglike = squeeze(cumsum(training_data.output_neglog_c{:}));
    training_loglike(L+1:end,:) = [];
    input_data{i} = training_loglike;
    nbData = [nbData, length(input_data{i})]; % so as to find out the max length
end
for n=1:length(input_data)
    align = spline(1:size(input_data{n},2), input_data{n}, linspace(1,size(input_data{n},2), max(nbData)));
    align_data = [align_data; align];%Resampling alignment
%     plot(align ,'LineWidth',1);
    hold on;
end

mean_loglike = mean(align_data);
if length(input_data)==1
    mean_loglike = align_data;
end
% plot(mean_loglike','b','LineWidth',3);
title('Caculate the mean log\_likelihood of training trials');
end