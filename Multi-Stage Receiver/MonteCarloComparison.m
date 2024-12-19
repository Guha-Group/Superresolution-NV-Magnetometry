%% MONTE-CARLO SIMULATIONS FOR BRIGHTNESS ESTIMATION

sigma = 1;
KAPPA = [0,0.1,0.2,0.3,0.4,0.5]; 
S = sigma*linspace(1e-2,1,100);
N = 1e5;

% run monte-carlo loop
num_mc_trials = 1e3;
KAPPA_EST = zeros(3,numel(KAPPA),numel(S),num_mc_trials); % container for estimates
for m1 = 1:numel(KAPPA)
    
    % assign parameter value
    k = KAPPA(m1);

    for m2 = 1:numel(S)
        % assign parameter value
        s = S(m2);

        for t = 1:num_mc_trials

            % Run Direct Imaging Experiment
            kappa_est_DD = DirectImagingExperiment(s,kappa,N,sigma);
            

            % Run Static Multi-Stage Receiver
            kappa_est_SR = StaticMultiStageReceiverExperiment(s,kappa,N,sigma);


            % Run Adaptive Multi-Stage Receiver
            kappa_est_AR = AdaptiveMultiStageReceiverExperiment(s,kappa,N,sigma);

            % assign results to the container
            KAPPA_EST(:,m1,m2,t) = [kappa_est_DD,kappa_est_SR,kappa_est_AR];
        end
    end
end


% Calculate MSE from monte-carlos
BIAS = sum(KAPPA_EST - KAPPA,4)/num_mc_trials;
MSE = sum((KAPPA_EST - KAPPA).^2,4)/num_mc_trials;

% Plot the MSE for each receiever
for r = 1:3
    for k = 1:numel(KAPPA)
        plot(MSE(r,:,:))
    end
end


