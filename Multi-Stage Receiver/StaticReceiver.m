% System Parameters
sigma = 1;          % diffraction limit
x0 = 0;             % geometric midpoint of the point sources
s = 0.3*sigma;      % half-separation of the point sources (must be greater than 0)
kappa = 0.3;        % relative brightness bias (must be between [-0.5,+0.5]

% Absolute Positions
x1 = x0-s;
x2 = x0+s;

% Photon Allocation
n1 = 1e5; % photons for first stage
n2 = 1e5; % photons for second stage
n3 = 1e5; % photons for third stage

%% STAGE 1: Direct Imaging

% Sources Balanced
b1 = 0.5;
b2 = 0.5;

% Simulate a measurement
X = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],n1,sigma);

% Estimate centroid and separation
[x0_est,s_est] = DirectImagingEstimates(X,sigma);

% Calculate moments for Gamma Prior
mu_est = (s_est/2/sigma)^2;
alpha = 8*n1*mu_est^2;
beta = 8*n1*mu_est;

%% STAGE 2: HG SPADE

% update pointing
x1 = x1-x0_est;
x2 = x2-x0_est;

% get HG mode sorting measurement
Q = SimulateHGSPADEMeasurement([x1,x2],[b1,b2],n2,sigma);

% update posterior distribution parameters
alpha = alpha + sum(Q,2);
beta = beta + n2;

%% STAGE 3: Direct Imaging

% set brightnesses
b1 = 0.5 - kappa;
b2 = 0.5 + kappa;

% simulate direct imaging measurement
X = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],n3,sigma);

% get MMSE for mu and s
mu_est = alpha/beta;
s_est = 2*sigma*sqrt(mu_est);

% estimate the brightness
kappa_est = mean(X,2)/2/s_est;
