%%% WARNING: THIS CODE GENERATES DATA FROM THE PRIOR, WHICH DOES NOT ALWAYS
%%% MEAN THAT THE REALIZATIONS HAVE DISTINGUISHABLE DYNAMICAL MODES. THIS
%%% CODE SIMPLY SERVES TO DEMONSTRATE HOW TO RUN THE HDPHMMDPinference
%%% SCRIPT.
% e.g runstuff(1, 1, 'F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox/saveDir')
  
%function runstuff(caseNumber,trial_vec,saveDir)
function run_pick_n_place()
clc;
close all;
caseNumber = 2;  %Observation mode type 
saveDir = 'F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox\saveDir';

%% Preprocess data 
ndemos = 1;
cellData = {}; 
data_struct = struct;
cellData = pick_n_place_file();
d = ndemos * size(cellData{1},2); % number of dimensitionality.

% use the 2nd order equation to smooth the original signal
 nbData = [];
 for i = 1: length(cellData) 
     cellData{i} = cellData{i}';
     nbData = [nbData, length(cellData{i})]; % so as to find out the max length
 end  
 sensor = zeros(d,max(nbData));
 initData = {};
 for n=1:length(cellData)
     initData = [initData; {cellData{i}}];
     cellData{n} = spline(1:size(cellData{n},2), cellData{n}, linspace(1,size(cellData{n},2), max(nbData))); %Resampling alignment
 end
sensor = [];
for i = 1: length(cellData)
    temp = cellData{i};
    sensor = [sensor; temp];
end 

%Smooth trajextory
radius = 50;
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
 
%Adjust each dim to mean 0
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
% generate the prior distribution by given model parameters, setting parameters and the observations (demonstrations) 
data_struct.obs = sensor;

data_struct.true_labels = ones(1, size(sensor, 2));

% --display-----  
fHandle = figure; 
for n = 1: ndemos
    subplot(3,1,1); plot(initData{n}');hold on;   %initial data
    title('Initial');  
    subplot(3,1,2); plot(cellData{n}');hold on;  %alignment
    title('Alignment');
end
subplot(3,1,3); plot(data_struct.obs'); hold on; plot(data_struct.true_labels'); 
title('Smooth | Transfer the mean = 0 and the covariance = Identity matrix');
hold off;
%}

%% Prior 
% obsModelType = {'Gaussian','AR','SLDS'}; % 3 kinds of observation model
switch caseNumber
    case 1
        obsModelType = 'Gaussian';
        priorType = 'IW-N';  %only the covariance is unknown;  prior on Gaussian N(mu_{k,j},Sigma_{k,j}) emissions (non-conjugate)
        sig0 = 10;  % covariance of the N(mu0,sig0) prior on mu_{k,j}
        meanSigma = eye(d);  % expected mean of IW(nu,nu_delta) prior on Sigma_{k,j}
        Kz = 10;  % truncation level of the DP prior on HMM transition distributions pi_k
        Ks = 10;  % truncation level of the DPMM on emission distributions pi_s
        
    case 2
        obsModelType = 'Gaussian';
        priorType = 'NIW';  %both the mean and the covariance are unknown: Normal-Inverse-Wishart distribution
                            %in NIW distribution the mean and the covariance are dependent on each other             
        kappa = 0.1;  % NIW(kappa,theta,delta,nu_delta) 
        %kappa:how strongly we believe the mean(affect the location) in this prior 
        %theta:the prior mean
        %delta:how strongly we believe the covariance(affect the shape) in this prior
        %nu_delta:propotional to the prior mean for covariance(meanSigma)
        
        meanSigma = eye(d); % the observation prior's covariance of gauessian distribution 
        
        %Kz: a number that exceeds the total number of expected HMM states
        Kz = 3; %truncation level of the DP(Dirichlet Process) prior on HMM transition distributions pi_k
        Ks = 3; %truncation level of the DPMM(Dirichlet Process Mixture Model) on emission distributions pi_s
        
    case 3
        obsModelType = 'AR';
        priorType = 'MNIW';  % prior on 
        r = 1;
        K = inv(diag([0.1*ones(1,d*r)]));
        meanSigma = eye(d);
        Kz = 5;
        
    case 4
        obsModelType = 'AR';
        priorType = 'MNIW-N';
        r = 2;
        K = inv(diag([0.1*ones(1,d*r)]));
        sig0 = 3;
        meanSigma = eye(d);
        Kz = 10;
        
    case 5
        obsModelType = 'AR';
        priorType = 'N-IW-N';
        r = 2;
        K = inv(diag([0.1*ones(1,d*r)]));
        sig0 = 3;
        meanSigma = eye(d);
        Kz = 10;
        
    case 6
        obsModelType = 'AR';
        priorType = 'ARD';
        r = 2;
        meanSigma = eye(d);        
        Kz = 10;
        
    case 7
        obsModelType = 'SLDS';
        dy = 2;
        priorType = 'MNIW';
        meanSigma = eye(d);
        K = inv(diag([0.1*ones(1,d)]));
        y_priorType = 'IW';
        y_var = 1;
        P0 = 1;
        Kz = 10;
        Kr = 10;
        
    case 8
        obsModelType = 'SLDS';
        dy = 1;
        priorType = 'MNIW';
        meanSigma = eye(d);
        K = inv(diag([0.1*ones(1,d)]));
        y_priorType = 'IW-N';
        y_var = 1;
        sig0_y = 1;
        P0 = 1;
        Kz = 10;
        Kr = 10;
        
    case 9
        obsModelType = 'SLDS';
        dy = 1;
        priorType = 'MNIW';
        meanSigma = eye(d);
        K = inv(diag([0.1*ones(1,d)]));
        y_priorType = 'NIW';
        y_var = 1;
        kappa_y = 1;
        P0 = 1;
        Kz = 10;
        Kr = 10;
        
    case 10
        obsModelType = 'SLDS';
        dy = 1;
        priorType = 'N-IW-N';
        K = inv(diag([0.1*ones(1,d)]));
        meanSigma = eye(d);
        sig0 = 10;
        y_priorType = 'IW';
        y_var = 1;
        P0 = 1;
        Kz = 10;
        
    case 11
        obsModelType = 'SLDS';
        dy = 1;
        priorType = 'N-IW-N';
        meanSigma = eye(d);
        K = inv(diag([0.1*ones(1,d)]));
        sig0 = 25;
        y_priorType = 'IW-N';
        y_var = 1;
        sig0_y = 25;
        P0 = 1;
        Kz = 10;
        Kr = 10;
        
end

%% Setting for inference:
clear settings
switch obsModelType
    case {'AR','SLDS'}
        settings.Kz = Kz;   % truncation level for mode transition distributions
        settings.Ks = 1;  % truncation level for mode transition distributions
        if strcmp(obsModelType,'SLDS') && ~strcmp(y_priorType,'IW')
            settings.Kr = Kr;  % truncation level for MoG measurement noise
        end
    case 'Gaussian'
        settings.Kz = Kz;   % truncation level for mode transition distributions
        settings.Ks = Ks;  % truncation level for mode transition distributions
end
settings.Niter = 2000;  % Number of iterations of the Gibbs sampler
settings.resample_kappa = 1;  % Whether or not to use sticky model
settings.seqSampleEvery = 100; % How often to run sequential z sampling
settings.saveEvery = 100;  % How often to save Gibbs sample stats
settings.storeEvery = 1;
settings.storeStateSeqEvery = 100;
settings.ploton = 1;  % Whether or not to plot the mode sequence while running sampler
settings.plotEvery = 1; %plot frequency
settings.plotpause = 0;  % Length of time to pause on the plot
settings.saveDir = saveDir;  % Directory to which to save files

%% Set Hyperparameters
clear model
% Type of dynamical system:
model.obsModel.type = obsModelType;
if strcmp(obsModelType,'AR')
    % Order of AR process:
    model.obsModel.r = r;
    m = d*r;
else
    m = d;
end

% Type of prior on dynamic parameters. Choices include matrix normal
% inverse Wishart on (A,Sigma) and normal on mu ('MNIW-N'), matrix normal
% inverse Wishart on (A,Sigma) with mean forced to 0 ('MNIW'), normal on A,
% inverse Wishart on Sigma, and normal on mu ('N-IW-N'), and fixed A,
% inverse Wishart on Sigma, and normal on mu ('Afixed-IW-N').  NOTE: right
% now, the 'N-IW-N' option is only coded for shared A!!!
model.obsModel.priorType = priorType;

switch model.obsModel.priorType
    case 'NIW'
        model.obsModel.params.M  = zeros([d 1]); %the prior mean
        model.obsModel.params.K =  kappa; % how strongly we believe in this prior
        
    case 'IW-N'
        % Mean and covariance for Gaussian prior on mean:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
    
    case 'MNIW'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);

        % Inverse covariance along rows of A (sampled Sigma acts as
        % covariance along columns):
        model.obsModel.params.K =  K(1:m,1:m);
        
    case 'MNIW-N'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);

        % Inverse covariance along rows of A (sampled Sigma acts as
        % covariance along columns):
        model.obsModel.params.K =  K(1:m,1:m);

        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));

    case 'N-IW-N'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]);
        model.obsModel.params.Lambda0_A = inv(kron(inv(K),meanSigma));

        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
        
    case 'Afixed-IW-N'
        % Set fixed A matrix:
        model.obsModel.params.A = A_shared;
        
        % Mean and covariance for mean of process noise:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
        
    case 'ARD'
        % Gamma hyperprior parameters for prior on precision parameter:
        model.obsModel.params.a_ARD = 10;
        model.obsModel.params.b_ARD = 0.01;
        
        % Placeholder for initializeStructs. Can I get rid of this?
        model.obsModel.params.M  = zeros([d m]);

        % Mean and covariance for mean of process noise:
        model.obsModel.params.zeroMean = 1;
end
        
% Degrees of freedom(model.obsModel.params.nu) and scale matrix(model.obsModel.params.nu_delta) for covariance of process noise:
% As model.obsModel.params.nu increases, the sampled matrices are more concentrated on the prior covariance matrix
model.obsModel.params.nu = 1000;%must be > d + 1: how strongly we believe the covariance(shape) in this prior
model.obsModel.params.nu_delta = (model.obsModel.params.nu-d-1)*meanSigma; % proportional to the prior mean for meanSigma (P20~2.12)ean for meanSigma (P20~2.12)
%model.obsModel.params.nu_delta = (model.obsModel.params.nu/(model.obsModel.params.nu-d-1))*meanSigma; % proportional to the prior m

if strcmp(obsModelType,'SLDS')
    % Degrees of freedom and scale matrix for covariance of measurement noise:
    model.obsModel.y_params.nu = 1000; %dy + 2;
    model.obsModel.y_params.nu_delta = (model.obsModel.y_params.nu-dy-1)*y_var*eye(dy);
    
    model.obsModel.y_priorType = y_priorType;
    
    switch model.obsModel.y_priorType
        case 'NIW'
            
            model.obsModel.y_params.M  = zeros([dy 1]);
            model.obsModel.y_params.K =  kappa_y;
            
        case 'IW-N'
            % Mean and covariance for Gaussian prior on mean:
            model.obsModel.y_params.mu0 = zeros(dy,1);
            model.obsModel.y_params.cholSigma0 = chol(sig0_y*eye(dy));
    end
    
    % Fixed measurement matrix:
    model.obsModel.params.C = [eye(dy) zeros(dy,d-dy)];
    
    % Initial state covariance:
    model.obsModel.params.P0 = P0*eye(d);
end

% Always using DP mixtures emissions, with single Gaussian forced by
% Ks=1...Need to fix.
model.obsModel.mixtureType = 'infinite';

% Sticky HDP-HMM parameter settings: weakly informative hyperprior, Gamma(1,0.01)prior on the concertration parapeters 
% alpha0_p_kappa0 = a_alpha / b_alpha;    % Gj concentration parameter
% gamma0 = a_gamma / b_gamma;    % G0 concentration parameter
model.HMMmodel.params.a_alpha=1;  % affects \pi_z
model.HMMmodel.params.b_alpha=0.01;
model.HMMmodel.params.a_gamma=1;  % global expected # of HMM states (affects \beta)
model.HMMmodel.params.b_gamma=0.01;
if settings.Ks>1
    model.HMMmodel.params.a_sigma = 1;
    model.HMMmodel.params.b_sigma = 0.01;
end
if isfield(settings,'Kr')
    if settings.Kr > 1
        model.HMMmodel.params.a_eta = 1;
        model.HMMmodel.params.b_eta = 0.01;
    end
end
model.HMMmodel.params.c=1000;  % self-transition proportion parameter. Beta(c, d) rho0 = c/(c+d)(In paper: rho = kappa/(kappa + alpha))
model.HMMmodel.params.d=1;
model.HMMmodel.type = 'HDP';

%%
%Generate data from the prior:
%{
  data_struct = generateData(model,settings,T);
  while length(unique(data_struct.true_labels))==1
      data_struct = generateData(model,settings,T);
  end
%}

% data_struct = PreprocessData(model,settings,length(sensor),sensor);
%  while length(unique(data_struct.true_labels))==1
%      data_struct = PreprocessData(model,settings,length(sensor),sensor);
%  end
 
%% -- -inference--------
trial_vec = [1:1];
for seq=1
    data_struct(1).test_cases = seq;
    for t=trial_vec
        settings.trial = t;  % Defines trial number, which is part of the filename used when saving stats
        [stateSeq, theta, hyperparams] =  HDPHMMDPinference(data_struct, model, settings);
        matfile = strcat('stateSeq',num2str(trial_vec),'.mat');
        save (matfile, 'stateSeq');
        %Calc the switchpoints and make lines on the traj graph
        switchpts = [];
        npts = length(stateSeq.z);
        last = stateSeq.z(1);
        for j=1:npts
            curr = stateSeq.z(j);
            if(curr ~= last)
                last = curr;
                switchpts = [switchpts, j];
            end      
        end
        figure(fHandle);
        gridxy(switchpts,'Linestyle','--')
   end
end
