%%% WARNING: THIS CODE GENERATES DATA FROM THE PRIOR, WHICH DOES NOT ALWAYS
%%% MEAN THAT THE REALIZATIONS HAVE DISTINGUISHABLE DYNAMICAL MODES. THIS
%%% CODE SIMPLY SERVES TO DEMONSTRATE HOW TO RUN THE HDPHMMDPinference
%%% SCRIPT.
% e.g runstuff(1, 1, 'F:\Matlab\Emily B. Fox\HDPHMM_HDPSLDS_toolbox/saveDir')
function train_state(trialID, training_state, training_dataset_path, LearnedModelPath)
clc;
close all;
global rootPath saveDir METHOD
global TRAINING_SIM_REAL TRAINING_SUCCESS_FAILURE TRAINING_PLOT_ON STATE 
global TRUNCATION_STATES TRUNCATION_COMPONENTS Niter
[data, R_State, foldname]= load_one_trial(trialID, training_dataset_path);
figure;
subplot(3,1,1);
plot(data');
grid MINOR
title('Original training signals')
[data, meanSigma] = data_preprocessing(data);
subplot(3,1,2);
plot(data');
grid MINOR
title('After data preprocessing');
gridxy(R_State(training_state:training_state + 1),'Linestyle','-','Color',[.5 .5 .5],'linewidth',2);

% specific step estimation
sensor = data(:,R_State(training_state):R_State(training_state + 1));

subplot(3,1,3);
plot(sensor');
grid MINOR
title(strcat('State: ', STATE(training_state)));

d = size(sensor,1); % number of dimensitionality.
%% Prior 
% obsModelType = {'Gaussian','AR','SLDS'}; % 3 kinds of observation model
switch METHOD
    case 'sHDPHMM'
        caseNumber = 2;  %Observation mode type 
    case 'HDPVARHMM(1)'
        caseNumber = 3;  %Observation mode type
    case 'HDPVARHMM(2)'
        caseNumber = 4;  %Observation mode type
    case 'HDPSLDS'
end
cd (saveDir); delete *; cd (rootPath);
switch caseNumber
    case 1
        obsModelType = 'Gaussian';
        priorType = 'IW-N';  %only the covariance is unknown;  prior on Gaussian N(mu_{k,j},Sigma_{k,j}) emissions (non-conjugate)
        sig0 = 10;  % covariance of the N(mu0,sig0) prior on mu_{k,j}
        meanSigma = eye(d);  % expected mean of IW(nu,nu_delta) prior on Sigma_{k,j}
        Kz = 4;  % truncation level of the DP prior on HMM transition distributions pi_k
        Ks = 4;  % truncation level of the DPMM on emission distributions pi_s       
    case 2
        obsModelType = 'Gaussian';
        priorType = 'NIW';  %both the mean and the covariance are unknown: Normal-Inverse-Wishart distribution
                            %in NIW distribution the mean and the covariance are dependent on each other             
        
        %kappa:how strongly we believe the mean(affect the location) in this prior 
        %mean:the prior mean
        %delta:how strongly we believe the covariance(affect the shape) in this prior
        %nu_delta:propotional to the prior mean for covariance(meanSigma)
        kappa = 0.1;  % NIW(kappa,mean,delta,nu_delta)    
        meanSigma = eye(d); % the observation prior's covariance of gauessian distribution 
        
        % for blocked sampler
        %Kz: a number that exceeds the total number of expected HMM states
        Kz = TRUNCATION_STATES;  %states : truncation level of the DP(Dirichlet Process) prior on HMM transition distributions pi_k
        Ks = TRUNCATION_COMPONENTS;   %components of mixture gaussian: truncation level of the DPMM(Dirichlet Process Mixture Model) on emission distributions pi_s
        
    case 3
        obsModelType = 'AR';
        priorType = 'MNIW';  % The conjugate matrix-normal inverse-Wishart prior for dynamic parameter () 
        r = 1;
        K = inv(diag([0.1*ones(1,d*r)]));
        meanSigma = eye(d);
        Kz = TRUNCATION_STATES;
        
    case 4
        obsModelType = 'AR';
        priorType = 'MNIW-N';
        r = 2;
        K = inv(diag([0.1*ones(1,d*r)]));
        sig0 = 1;  %default: 3
        meanSigma = eye(d);
        Kz = TRUNCATION_STATES;
        
    case 5
        obsModelType = 'AR';
        priorType = 'N-IW-N';
        r = 2;
        K = inv(diag([0.1*ones(1,d*r)]));
        sig0 = 3;
        meanSigma = eye(d);
        Kz = TRUNCATION_STATES;
        
    case 6
        obsModelType = 'AR';
        priorType = 'ARD';
        r = 2;
        meanSigma = eye(d);        
        Kz = TRUNCATION_STATES;
        
    case 7
        obsModelType = 'SLDS'; %switching linear dynamical system
        dy = 2;
        priorType = 'MNIW';
        meanSigma = eye(d);
        K = inv(diag([0.1*ones(1,d)]));
        y_priorType = 'IW';
        y_var = 1;
        P0 = 1;
        Kz = TRUNCATION_STATES;
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

% generate the prior distribution by given model parameters, setting parameters and the observations (demonstrations) 
data_struct.obs = sensor;
data_struct.true_labels =  ones(1, size(sensor, 2));

%% Setting for inference:
clear settings
switch obsModelType
    case {'AR','SLDS'}
        settings.Kz = Kz;   % truncation level for mode transition distributions
        settings.Ks = 1;    % truncation level for mode transition distributions
        if strcmp(obsModelType,'SLDS') && ~strcmp(y_priorType,'IW')
            settings.Kr = Kr;           % truncation level for MoG measurement noise
        end
    case 'Gaussian'
        settings.Kz = Kz;               % truncation level for mode transition distributions
        settings.Ks = Ks;               % truncation level for mode transition distributions
end
settings.Niter              = Niter;     % Number of iterations of the Gibbs sampler
settings.resample_kappa     = 1;        % Whether or not to use sticky model
settings.seqSampleEvery     = 1;       % How often to run sequential z sampling
settings.saveEvery          = 1000;       % How often to save Gibbs sample stats
settings.storeEvery         = 1000;
settings.storeStateSeqEvery = 1000;
settings.ploton             = TRAINING_PLOT_ON;        % Whether or not to plot the mode sequence while running sampler
settings.plotEvery          = 1;        %plot frequency
settings.plotpause          = 0;        % Length of time to pause on the plot
settings.saveDir            = saveDir;  % Directory to which to save files

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
% inverse Wishart on (A,Sigma) and normal on mu ('MNIW-N'),
% matrix normal inverse Wishart on (A,Sigma) with mean forced to 0 ('MNIW'), 
% normal on A,inverse Wishart on Sigma, and normal on mu ('N-IW-N'), and 
% fixed A, inverse Wishart on Sigma, and normal on mu ('Afixed-IW-N').  
% NOTE: right now, the 'N-IW-N' option is only coded for shared A!!!
model.obsModel.priorType = priorType;

switch model.obsModel.priorType
    case 'NIW'
        model.obsModel.params.M  = zeros([d 1]);    % the prior mean
        model.obsModel.params.K =  kappa;           % how strongly we believe the mean in this prior
        
    case 'IW-N'
        % Mean and covariance for Gaussian prior on mean:
        model.obsModel.params.mu0 = zeros(d,1);
        model.obsModel.params.cholSigma0 = chol(sig0*eye(d));
    
    case 'MNIW'
        % Mean and covariance for A matrix:
        model.obsModel.params.M  = zeros([d m]); %mean

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
switch obsModelType
    case 'Gaussian'
        model.obsModel.params.nu = 1000;                                            %must be > d + 1: how strongly we believe the covariance(shape) in this prior
   case 'AR'
        model.obsModel.params.nu = d + 2;                                            %must be > d + 1: how strongly we believe the covariance(shape) in this prior
   case 'SLDS'
        model.obsModel.params.nu = d + 2;                                            %must be > d + 1: how strongly we believe the covariance(shape) in this prior
end
model.obsModel.params.nu_delta = (model.obsModel.params.nu-d-1)*meanSigma;  % proportional to the prior mean for meanSigma (P20~2.12)ean for meanSigma (P20~2.12)

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
                                                  % Sticky HDP-HMM hyperparameters settings: 
                                                  % Gamma(1,0.01)prior on the (alpha0 + kappa0)
                                                  % alpha0_p_kappa0 = a_alpha / b_alpha; 
model.HMMmodel.params.a_alpha=1;                  % default: 1 affects \pi_z
model.HMMmodel.params.b_alpha=0.01;               % default:0.01

                                                  % Gamma(1,0.01)prior on the concertration parapeters :gamma0 = a_gamma / b_gamma;    % G0 concentration parameter
model.HMMmodel.params.a_gamma=1;                  % default: 1 global expected # of HMM states (affects \beta)
model.HMMmodel.params.b_gamma=0.01;               % default:0.01
if settings.Ks>1
    model.HMMmodel.params.a_sigma = 1;            % default: 1
    model.HMMmodel.params.b_sigma = 0.01;         % default: 0.01
end
if isfield(settings,'Kr')
    if settings.Kr > 1
        model.HMMmodel.params.a_eta = 1;
        model.HMMmodel.params.b_eta = 0.01;
    end
end

                                                  % self-transition proportion parameter rho modeled by a Beta
switch obsModelType
    case 'Gaussian'
        model.HMMmodel.params.c=100;              %Beta(c, d) rho0 = c/(c+d)(In paper: rho = kappa/(kappa + alpha))
    case 'AR'
        model.HMMmodel.params.c=10;               %Beta(c, d) rho0 = c/(c+d)(In paper: rho = kappa/(kappa + alpha))
    case 'SLDS'
        model.HMMmodel.params.c=10;               %Beta(c, d) rho0 = c/(c+d)(In paper: rho = kappa/(kappa + alpha))
end
model.HMMmodel.params.d=1;                        %default: 1
model.HMMmodel.type = 'HDP';

%% -- -inference--------
settings.trial = 1;  % Defines trial number, which is part of the filename used when saving stats
[stateSeq,                ...
    theta,                ...
    hyperparams,          ...
    neglog_c,             ... 
    total_log_likelihood, ...
    loglike_normalizer,   ... 
    dist_struct,obsModel, ...
    data_struct] =  HDPHMMDPinference(data_struct, model, settings,char(STATE(training_state)),trialID);  
%% Postprocessing
cd (LearnedModelPath);
LearnedModelFile = strcat(TRAINING_SIM_REAL, '_', TRAINING_SUCCESS_FAILURE, '_',char(STATE(training_state)),'_',foldname,'.mat');

save (LearnedModelFile, 'stateSeq', ...
                'theta', ...
                'hyperparams', ...
                'neglog_c', ...
                'total_log_likelihood',...
                'loglike_normalizer',...
                'dist_struct',...
                'obsModel',...
                'data_struct');
disp('................................................................................................................................');
