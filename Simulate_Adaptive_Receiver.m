function [theta_mle, pho_tot, alpha_opt, schedule_opt] =...
    Simulate_Adaptive_Receiver(theta,adapt_schedule,IMG,ODMR)
    % theta=[s,d1,d2]   : Ground truth parameter values
    % adapt_schedule    : True/False flag for also performing microwave scheduling optimization
    % IMG               : image system struct
    % ODMR              : ODMR protocol struct
    
    % unpack parameters
    s  = theta(1);
    d1 = theta(2);
    d2 = theta(3);
    
    % setup imaging system
    sigma = IMG.sigma;
    integration_time = IMG.integration_time;
    q = IMG.q;
    x = IMG.x;
    
    % setup ODMR system
    eta0 = ODMR.eta0;
    chi = ODMR.chi;
    w0 = ODMR.w0;
    schedule = ODMR.schedule;
    linewidth = ODMR.linewidth;
    omega = ODMR.omega;

    % temporal time-step
%   dt = integration_time/numel(ODMR.omega); % discrete time steps
    dt = integration_time/50; % discrete time steps
    t = 0;                      % initialize clock

    %% STAGE 1: SPADE
    P1 = 0.5*(HG_prob(-s/2,sigma,q) + HG_prob(+s/2,sigma,q));   % probability of detection
    Q = zeros(size(q));                                         % container for stage 1 measurement
    alpha_t = 1;                                                % adaptive estimator of optimal switching parameter
    eta1 = 2*eta0;
    
    % adaptive switching
    while (t < integration_time)  &&  (t/integration_time < alpha_t)
        
        % collect measurement at time interval
        Q_t = SimulateMeasurement(P1, eta1*dt);

        if any(Q_t>0)
            
            % add measurement results to cumulative stage 1 measurements
            Q = Q + Q_t;

            % estimate the separation
            s_hat = SPADE_MLE(Q,q,sigma);

            % add a little bias to the estimator to avoid div by 0s
            if s_hat == 0
                s_hat = 1e-6*rand(1)*IMG.sigma;
            end

            % optimze cost function to get alpha_t
            alpha_t = AdaptiveSwitchingOptimization(IMG, ODMR,s_hat);

        end
        
        % increment time
        t = t + dt;

    end
    
    % adaptive switching time
    t1 = t;
    alpha_opt = t1/integration_time;

    %% STAGE 2: ODMR
    [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);
    eta2 = I1 + I2;
    b1 = I1./eta2;
    b2 = I2./eta2;

    % Get direct imaging ODMR probability distribution 
    p1 = DD_prob(-s/2,sigma,x); p1 = p1/sum(p1);
    p2 = DD_prob(+s/2,sigma,x); p2 = p2/sum(p2);
    P2 = b1.*p1 + b2.*p2;

    % adaptive scheduling
    if ~adapt_schedule
        % collect measurement at time interval
        t2 = integration_time-t1;
        X = SimulateMeasurement(P2, t2.*eta2.*schedule);
        schedule_opt = schedule;
    else
        
        X = zeros([numel(omega),numel(x)]); % container for direct imaging measurements
        omega_id_t = 1;                     % adaptive drive frequency index
        omega_id = [];
    
        while (t < integration_time)
    
            % collect measurement at time interval
            X_t = SimulateMeasurement(P2(omega_id_t,:), eta2(omega_id_t)*dt);
    
            if any(X_t>0,"all")
    
                % add measurement to cumulative stage 2 measurements
                X(omega_id_t,:) = X(omega_id_t,:) + X_t;
    
                % estimate parameters
                theta_hat = TWOSTAGE_MLE(Q,q,X,x,IMG,ODMR);
    
                % optimize the drive frequency
                omega_id_t = AdaptiveSchedulingOptimization(IMG,ODMR,theta_hat);
    
            end
            % add current frequency index to the cumulative schedule 
            omega_id = [omega_id, omega_id_t];
    
            % increment time
            t = t + dt;
        end
        
        % recover the optimal schedule
        [gc,gr]=groupcounts(omega_id');
        schedule_opt = zeros(size(omega));
        schedule_opt(gr) = gc/sum(gc);
    end


    % final estimation and photon count
    theta_mle = TWOSTAGE_MLE(Q,q,X,x,IMG,ODMR);
    pho_tot = sum(Q,'all') + sum(X,'all');    
end