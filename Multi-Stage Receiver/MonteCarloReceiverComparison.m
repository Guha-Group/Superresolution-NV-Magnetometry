%% MONTE-CARLO SIMULATIONS FOR BRIGHTNESS ESTIMATION
sigma = 1;
x0 = 0;
S = sigma*linspace(1e-2,1,100);
KAPPA = [0,0.1,0.2,0.3,0.4,0.5]; 
N = 1e5;

% special arguments for each receiver

%% DIRECT IMAGING RECEIVER 
max_iterations = 100;

%% STATIC RECEIVER
alpha_1 = 1/3; 
alpha_2 = 1/3; 

%% ADAPTIVE RECEIVER
pho_step = 100;

% run monte-carlo loop
num_mc_trials = 1e4;
PARAM_EST = zeros(3,2,numel(KAPPA),numel(S),num_mc_trials);

parpool(96)
for m1 = 1:numel(KAPPA)
    
    % assign parameter value
    kappa = KAPPA(m1);

    for m2 = 1:numel(S)
        
        % assign parameter value
        s = S(m2);

        parfor t = 1:num_mc_trials

            % Run Direct Imaging Experiment
            param_est_DD = DirectImagingReceiver(x0,s,kappa,N,sigma,max_iterations);

            % Run Static Multi-Stage Receiver
            param_est_SR = StaticMultiStageReceiver(x0,s,kappa,N,sigma,alpha_1,alpha_2);

            % Run Adaptive Multi-Stage Receiver
            %param_est_AR = AdaptiveMultiStageReceiver(x0,s,kappa,N,sigma,pho_step);

            % assign results to the container
            % PARAM_EST(:,:,m1,m2,t) = [param_est_DD,param_est_SR,param_est_AR];
            PARAM_EST(:,:,m1,m2,t) = [param_est_DD,param_est_SR];

        end
    end
end


save('ReceiverComparison_DD_SR.mat','PARAM_EST')


% Calculate MSE from monte-carlos
BIAS = sum(KAPPA_EST - KAPPA,4)/num_mc_trials;
MSE = sum((KAPPA_EST - KAPPA).^2,4)/num_mc_trials;


