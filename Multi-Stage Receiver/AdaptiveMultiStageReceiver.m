function params_out = AdaptiveMultiStageReceiver(x0,s,kappa,N,sigma,pho_step)

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
    
    addpath('utils/')
 
    % Scene Parameters
    % x0 - geometric midpoint of the point sources
    % s  - half-separation of the point sources (must be greater than 0)
    % kappa - relative brightness bias (must be between [-0.5,+0.5]

    % System Parameters
    % sigma     - diffraction limit
    % N         - photon budget
    % pho_step  - number of photons detected betweeen adaptive updates

    
    %% ADAPTIVE ESTIMATION:
        
    % Absolute Positions
    x1 = x0-s;
    x2 = x0+s;
    
    % Measurement Parameters
    n1 = 0;         % number of photons allocated to stage 1 (as of current time)
    n2 = 0;         % number of photons allocated to stage 2 (as of current time)
    
    % load in optimal switch data for first stage
    load('SwitchingLookupTables.mat','N_range','s_range','opt_switch_frac')
    mu_range = (s_range/2/sigma).^2;
    min_mu = min(mu_range); max_mu = max(mu_range);
    
    
    %% STAGE 1: DIRECT IMAGING
    
    % Balance the source brightnesses
    b1 = 0.5;   % brightness of source 1
    b2 = 0.5;   % brightness of source 2
    
    % Initialize adaptive paremeters
    mu_est = 0; % rate parameter mu = (s/2/sigma^2)
    n1_opt = N; % optimal photon allocation the first stage
    X = [];     % container for direct imaging photon arrival positions
    
%   while (n1<=N) && (n1 < n1_opt || mu_est == 0 ) 
    while (n1<=N) && (n1 < n1_opt || n1 < 1e2 * (1+4*mu_est)/(4*mu_est))   
        
        % simulate direct imaging measurement
        X = [X,SimulateDirectImagingMeasurement([x1,x2],[b1,b2],pho_step,sigma)];
    
        % add the new photons to the stage 1 counter
        n1 = n1 + pho_step;
    
        % estimate midpoint x0, half separation s, and rate parameter mu
        %[x0_est,s_est] = DirectImagingMLEs(X,sigma);
        
        % EM estimation of source positions
        position_flag = 1;
        brightness_flag = 0;
        [x1_est,x2_est,~,~] = ExpectationMaximizationDD(X,sigma,100,position_flag,[nan,nan],brightness_flag,[0.5,0.5]);
        x0_est = (x2_est + x1_est)/2;
        s_est  = (x2_est - x1_est)/2;
        
        mu_est = (s_est/2/sigma)^2;
        
        % If the rate parameter (separation) is greater than zero...
        if (mu_est > 0)
    
            % calculate optimal stage 3 photons
            n3_opt = max(0,round(N-sqrt(N/mu_est)));
            
            % calculate photons remaining for the first and second stages
            N12 = N-n3_opt;
    
            % optimize switching time between the first and second stage
            % (Grace's Taylor Expansion Approach)
            %[n1_opt,~] = OptimizeFirstSwitching(s_est,N12,N12/10,sigma);
            
            [~,nearest_N_id] = min(abs(N_range - N12));
            n1_opt = interp1(mu_range,N12*opt_switch_frac(nearest_N_id,:),mu_est,'linear');
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
    %{
    temp_mu = linspace(0,alpha/beta^2,1000);
    plot(temp_mu,gampdf(temp_mu,alpha,beta));
    %}
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
        
        %{
        hold on
        temp_mu = linspace(0,alpha_post/beta_post^2,1000);
        plot(temp_mu,gampdf(temp_mu,alpha_post,beta_post));
        hold off
        %}

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
    %kappa_est = mean(XX)/(2*s_est);

    % Run EM with positions constrained
    brightness_flag = 1;
    position_flag = 0;
    [~,~,b1_est,b2_est] = ExpectationMaximizationDD(XX,sigma, 100, position_flag, s_est*[-1,1], brightness_flag, [nan,nan]);
    kappa_est = (b2_est-b1_est)/2;

    % collect parameter estimates
    params_out = [x0_est; s_est; kappa_est];
end

