function theta_mle = TWOSTAGE_MLE(Q, q, X, x,IMG,ODMR)
    % Q     : [1,N] measurement vector containing the number of photon arrivals in each mode of the SPADE measurement
    % q     : [1,N] vector of mode indices
    % X     : [W,M] measurement matrix containing number of photon arrival at each focal plane position for each odmr frequency
    % x     : [1,M] the direct imaging detector positions

    % ODMR and system variables
    sigma = IMG.sigma;
    eta0 = ODMR.eta0;
    chi = ODMR.chi;
    omega = ODMR.omega;
    w0 = ODMR.w0;
    linewidth = ODMR.linewidth;
    

    % ADAM optimizer gradient ascent parameters
    num_inits = IMG.num_inits;          % number of times to randomly initalize the gradient ascent algorithm (robustness to local maxima)
    num_steps = IMG.num_steps;          % number of gradient ascent iterations
    learning_rate = IMG.learning_rate;  % step size in the gradient.
    beta_1 = IMG.beta_1;                % gradient momentum param 1
    beta_2 = IMG.beta_2;                % gradient momentum param 2
    epsilon = IMG.epsilon;              % a small constant
    
    % containers for log-likelihood and intermediate estimates
    mle = zeros([3,num_inits,num_steps+1]); 
    loglike = zeros(num_inits,num_steps+1);
    
    for i = 1:num_inits

        % ADAM optimizer mean gradient track
        mt = zeros(3,1);
        vt = zeros(3,1);

        % (pseudo)-random initialization for MLE gradient ascent 
        [d1,d2] = DetuningPriorRnd(ODMR,1);
        theta_mle = [SPADE_MLE(Q,q,sigma), d1,d2]';
        mle(:,i,1) = theta_mle; 

        for j = 1:num_steps
            
            % expand parameters
            s = theta_mle(1);
            d1 = theta_mle(2);
            d2 = theta_mle(3);
            
            %% STAGE 1 GRADIENT:

            % probability over sample space
            P1 = 0.5*(HG_prob(-s/2,sigma,q) + HG_prob(+s/2,sigma,q));
            
            % partial derivative wrt to s
            dLogP1_ds = .5 * .5 * sum(Q .* (-d_HG_prob(-s/2,sigma,q) + d_HG_prob(+s/2,sigma,q)) ./ P1 ,2);

            % gradient
            grad1 = [dLogP1_ds,0,0]';

            %% STAGE 2 GRADIENT:
            % compute I1, I2
            [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);
            
            % compute photon rate over sample space
            PI = (I1 .* DD_prob(-s/2,sigma,x) + I2 .* DD_prob(+s/2,sigma,x));

            % probability over sample space
            P2 = PI./(I1+I2);

            % compute derivatives of I1 (I2) wrt to detuning d1 (d2)
            [dI1,dI2] = d_ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);
            
            % partial derivative wrt s
            dLogP_ds = .5 * sum(X.*( -I1 .* d_DD_prob(-s/2,sigma,x) + I2 .* d_DD_prob(+s/2,sigma,x) ) ./ PI, [1,2]);
            
            % partial derivative wrt d1
            dLogP_dd1 = sum( X.*(dI1 .* DD_prob(-s/2,sigma,x) ./ PI - dI1./(I1+I2)), [1,2]);

            % partial derivative wrt d2
            dLogP_dd2 = sum( X.*(dI2 .* DD_prob(+s/2,sigma,x) ./ PI - dI2./(I1+I2)), [1,2]);

            % gradient
            grad2 = [dLogP_ds, dLogP_dd1, dLogP_dd2]';
            

            %% LOG-LIKELIHOOD AT CURRENT ESTIMATE 
            loglike(i,j) = sum(Q.*log(P1),2) + sum(X.*log(P2),[1,2]);

            %% PERFORM GRADIENT ASCENT UPDATE STEP USING ADAM OPTIMIZER
            grad = grad1 + grad2;
            
            % ADAM optimizer gradient updates
            mt = beta_1*mt + (1-beta_1)*grad;
            vt = beta_2*vt + (1-beta_2)*grad.^2;
            mt_bar = mt/(1-beta_1);
            vt_bar = vt/(1-beta_2);
            
            % parameter update step
            theta_mle = theta_mle + (learning_rate./(sqrt(vt_bar)+epsilon)).* mt_bar; 

            % enforce positivity constraint on parameters
            theta_mle = abs(theta_mle);

            % store the MLE for the current random initializaion
            mle(:,i,j+1) = theta_mle;
        end

        % expand parameters
        s = theta_mle(1);
        d1 = theta_mle(2);
        d2 = theta_mle(3);
        
        % compute first stage probability
        P1 = 0.5*(HG_prob(-s/2,sigma,q) + HG_prob(+s/2,sigma,q));

        % compute second stage probabtility
        [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);
        P2 = (I1 .* DD_prob(-s/2,sigma,x) + I2 .* DD_prob(+s/2,sigma,x)) ./ (I1+I2);     

        % store the log likelihood of the current MLE random initialization 
        loglike(i,j+1) = sum(Q.*log(P1),2) + sum(X.*log(P2),[1,2]);

    end

    % choose MLE candidate with greatest log likelihood
    [~,k] = max(loglike,[],'all');
    [ki,kj] = ind2sub(size(loglike),k);
    theta_mle = mle(:,ki,kj);
end


%{
function theta_mle = TWOSTAGE_MLE(Q, q, X, x,IMG,ODMR)
    
    % unpack variables
    omega = ODMR.omega;
    w0 = ODMR.w0;
    linewidth =ODMR.linewidth;
    eta0 = ODMR.eta0;
    chi = ODMR.chi;
    sigma=IMG.sigma;

    % estimate separation from SPADE measurements
    if any(Q>0)
        s_est = SPADE_MLE(Q,q,sigma);
    else
        s_est = 1e-6*sigma;
    end

    % filter ODMR frequencies with no photon arrivals
    omega_nz = sum(X,2)>0;
    omega = omega(omega_nz);
    X = X(omega_nz,:);

    % get predicted proibability distributions
    p1 = DD_prob(+s_est/2,sigma,x);
    p2 = DD_prob(-s_est/2,sigma,x);

    % fit the brightnesses
    EM_iters = 1e4; % maximum number of EM iterations
    b_est = .5* ones(size(X,1),1);
    for k = 1:EM_iters
        b_est = sum(X.* (b_est.*p1)./(b_est.*p1+ (1-b_est).*p2),2)./sum(X,2);
        b_est(isnan(b_est))=0;
    end
    % estimates of emission spectra
    I1_est = b_est.*sum(X,2);
    I2_est = (1-b_est).*sum(X,2);

    % Define Lorentzian and its derivative wrt to the detuning
    Lorentzian = @(omega, d) 1 ./ (1 + (omega - w0 - d).^2 / linewidth^2);
    
    % initialize guesses of d1, d2
    [d1_est,d2_est] = DetuningPriorRnd(ODMR,1);
    
    % fit model to brightnesses
    model = @(d,freq) eta0 * (1 - chi/2 * (Lorentzian(freq,-d) + Lorentzian(freq,+d)));
    ft = fittype(model,'indep','freq');
    s1_params = fit(omega,I1_est,ft,'lower',min(omega)-w0,'upper',max(omega)-w0,'StartPoint',d1_est);
    s2_params = fit(omega,I2_est,ft,'lower',min(omega)-w0,'upper',max(omega)-w0,'StartPoint',d2_est);

    d1_est = s1_params.d;
    d2_est = s2_params.d;

    theta_mle = [s_est,d1_est,d2_est]';
end
%}
