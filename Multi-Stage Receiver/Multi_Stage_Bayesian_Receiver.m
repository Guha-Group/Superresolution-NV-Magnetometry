%% MULTI-STAGE RECEIEVER FOR 2-NV ESTIMATION
% Description:
% Simulates a three stage receiver for estimating the brightness of two
% point nitrogen-vacancy emitters that reside below the diffraction limit
% of resolution for an imaging system with a gaussian PSF. The stage
% sequence is:
% STAGE 1: DIRECT IMAGING (balanced brightnesses)
% - estimates the geometric midpoint between the sources
% - pre-estimates the separation between the sources
%
% SECOND STAGE: HERMITE-GAUSS SPADE (balanced brightnesses)
% - refines estimate of the separation between the sources
%
% THIRD STAGE: DIRECT IMAGING (unbalanced brightnesses)
% - estimates the relative brightness of the sources
% 
% 
% Author: Nico Deshler
% Date: 12/12/2024

clear all

%% ADAPTIVE ESTIMATION:

% System Parameters
sigma = 1;          % diffraction limit
x0 = 0;             % geometric midpoint of the point sources
s = 0.3*sigma;      % half-separation of the point sources (must be greater than 0)
kappa = 0.3;        % relative brightness bias (must be between [-0.5,+0.5]
mu = (s/2/sigma)^2;

% Absolute Positions
x1 = x0-s;
x2 = x0+s;

% Measurement Parameters
N = 1e5;        % total photon budget
pho_step = 1e2; % how many photons to detect betweeen adaptive updates
n1 = 0;         % number of photons allocated to stage 1 (as of current time)
n2 = 0;         % number of photons allocated to stage 2 (as of current time)

% load in optimal switch data for first stage
%{
load('SwitchingLookup.mat')
mu_range = (s_range/2/sigma).^2;
min_mu = min(mu_range); max_mu = max(mu_range);
%}

%% STAGE 1: DIRECT IMAGING

% Balance the source brightnesses
b1 = 0.5;   % brightness of source 1
b2 = 0.5;   % brightness of source 2

% Initialize adaptive paremeters
mu_est = 0; % rate parameter mu = (s/2/sigma^2)
n1_opt = N; % optimal photon allocation the first stage
X = [];     % container for direct imaging photon arrival positions

while (mu_est == 0 || n1_opt > n1) && (n1<=N)
    
    % simulate direct imaging measurement
    X = [X,SimulateDirectImagingMeasurement([x1,x2],[b1,b2],pho_step,sigma)];

    % add the new photons to the stage 1 counter
    n1 = n1 + pho_step;

    % estimate midpoint x0, half separation s, and rate parameter mu
    [x0_est,s_est] = DirectImagingEstimates(X,sigma);
    mu_est = (s_est/2/sigma)^2;
    
    % If the rate parameter (separation) is greater than zero...
    if (mu_est > 0)

        % calculate optimal stage 3 photons
        n3_opt = max(0,round(N-sqrt(N/mu_est)));
        
        % calculate photons remaining for the first and second stages
        N12 = N-n3_opt;

        % optimize switching time between the first and second stage
        % (Grace's Taylor Expansion Approach)
        [n1_opt,~] = OptimizeFirstSwitching(s_est,N12,sigma);
        
        %n1_opt = interp1(mu_range,N12*opt_switching_fraction,mu_est,'linear');
    else
        n1_opt = N;
    end
end

% update pointing
x1 = x1 - x0_est;
x2 = x2 - x0_est;

% set Gamma prior hyperparameters for rate parameter mu
alpha = n1 * mu_est^2;
beta = n1 * mu_est;
temp_mu = linspace(0,alpha/beta^2,1000);
plot(temp_mu,gampdf(temp_mu,alpha,beta));

%% STAGE 2: SPADE
n2_opt = N12 - n1_opt;
Q = [];     % container for SPADE measurement outcomes (photon indices)
while (N-n3_opt) > (n1+n2) && n2_opt > n2
    
    % simulate SPADE measurement
    Q = [Q,SimulateHGSPADEMeasurement([x1,x2],[b1,b2],pho_step,sigma)];

    % add the new photons to the stage 2 counter
    n2 = n2 + pho_step;

    % update hyperparameters for gamma prior on rate mu 
    alpha_post = alpha + sum(Q);
    beta_post = beta + numel(Q);
    
    hold on
    temp_mu = linspace(0,alpha_post/beta_post^2,1000);
    plot(temp_mu,gampdf(temp_mu,alpha_post,beta_post));
    hold off

    % estimate mu
    mu_est = alpha_post/beta_post;

    % evaluate optimal switching time for third stage
    n3_opt = N-sqrt(N/mu_est);

    % evaluate optimal switching time for second stage given first stage
    % counts and optimal third stage switching
    n2_opt = N - n3_opt - n1;
end


%% STAGE 3: Direct Imaging (Brightness Estimation Step)

% set brightnesses
b1 = 0.5 - kappa;
b2 = 0.5 + kappa;

% simulate direct imaging measurement
n3 = N - (n1 + n2);
XX = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],n3,sigma);

% estimate brightness detuning
s_est = 2*sigma*sqrt(mu_est);
kappa_est = mean(XX)/(2*s_est);

%% FUNCTIONS
function [x0_est, s_est] = DirectImagingEstimates(X,sigma)
    
    % ML estimator equations
    x0_est = mean(X);
    s_ML_transcendental = @(s_ML, x0) s_ML - mean( tanh((X-x0).*s_ML./sigma^2).*(X-x0),2);
    x0_ML_transcendental = @(x0_ML, s) x0_ML - mean(X + s * tanh((X-x0_ML) .* s / (sigma^2)),2);
    
    % alternate estimators twice to get convergence
    for k = 1:2
        Ls = @(s_ML) s_ML_transcendental(s_ML, x0_est);
        s_est = abs(fzero(Ls,sigma));
        Lx = @(x0_ML) x0_ML_transcendental(x0_ML,s_est);
        x0_est = fzero(Lx,0);
    end
end

function Q = SimulateHGSPADEMeasurement(xs,bs,num_pho,sigma)
    num_sources = numel(bs);    % total number of sources
    ns = round(bs*num_pho);     % number of photons from each source

    Q = zeros(1,num_pho);       % container for photon mode indices
    kk=1;                       % index counter
    for k = 1:num_sources
        lambda = (xs(k)/2/sigma)^2;
        Q(kk:(kk-1+ns(k))) = poissrnd(lambda,1,ns(k));
        kk = kk + ns(k);
    end
end

function X = SimulateDirectImagingMeasurement(xs,bs,num_pho,sigma)
    % Description: Returns the arrival location of photons on an ideal
    % photo-detector.
    num_sources = numel(bs);    % total number of sources
    ns = round(bs*num_pho);     % number of photons from each source
    
    X = zeros(1,num_pho);         % container for photon arrival locations
    kk=1;                       % index counter
    for k = 1:num_sources
        X(kk:(kk-1+ns(k))) = sigma*randn(1,ns(k))+xs(k);
        kk = kk+ns(k);
    end
end

