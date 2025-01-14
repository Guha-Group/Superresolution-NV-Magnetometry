function params_out = StaticMultiStageReceiver(x0,s,kappa,N,sigma,alpha_1,alpha_2)

    addpath('utils/')
 
    % Scene Parameters
    % x0 - geometric midpoint of the point sources
    % s  - half-separation of the point sources (must be greater than 0)
    % kappa - relative brightness bias (must be between [-0.5,+0.5]

    % System Parameters
    % sigma     - diffraction limit
    % N         - photon budget
    % alpha_1   - fraction of photon budget allocated to first stage
    % alpha_2   - fraction of photon budget allocated to second stage
    
    % Absolute Positions
    x1 = x0-s;
    x2 = x0+s;
    
    % Photon Allocation
    n1 = round(N*alpha_1); % photons for first stage
    n2 = round(N*alpha_2); % photons for second stage
    n3 = N-n1-n2;          % photons for third stage
    
    %% STAGE 1: Direct Imaging
    
    % Sources Balanced
    b1 = 0.5;
    b2 = 0.5;
    
    % Simulate a measurement
    X = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],n1,sigma);
    
    % Estimate centroid and separation
    %[x0_est,s_est] = DirectImagingMLEs(X,sigma);
    position_flag = 1;
    brightness_flag = 0;
    [x1_est,x2_est,~,~] = ExpectationMaximizationDD(X,sigma, 100, position_flag, [nan,nan], brightness_flag, [0.5,0.5]);
    x0_est = (x2_est + x1_est)/2;
    s_est =  (x2_est - x1_est)/2;

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
    % kappa_est = mean(X,2)/2/s_est;
    position_flag = 0;
    brightness_flag = 1;
    [~,~,b1_est,b2_est] = ExpectationMaximizationDD(X,sigma, 100, position_flag, s_est*[-1,1], brightness_flag, [nan,nan]);
    kappa_est = (b2_est - b1_est)/2;

    % collect parameter estimates
    params_out = [x0_est; s_est; kappa_est];
end
