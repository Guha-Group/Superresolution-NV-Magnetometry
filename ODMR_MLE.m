function theta_mle = ODMR_MLE(X, x, IMG, ODMR)
    % compute MLE from the simulated ODMR measurement via gradient ascent
    % on the likelihood
    
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
        
        % random initialization of the parameters
        [d1,d2] = DetuningPriorRnd(ODMR,1);
        theta_mle = [abs(normrnd(0,sigma/3)), d1,d2]';
        mle(:,i,1) = theta_mle; 

        for j = 1:num_steps
                
                % expand parameters
                s = theta_mle(1);
                d1 = theta_mle(2);
                d2 = theta_mle(3);
                
                % compute I1, I2
                [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);

                % compute derivatives of I1 (I2) wrt to detuning d1 (d2)
                [dI1,dI2] = d_ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);
                
                % compute probability over sample space
                PI = (I1 .* DD_prob(-s/2,sigma,x) + I2 .* DD_prob(+s/2,sigma,x));
                P = PI./(I1+I2);
                
                % partial derivative wrt s
                dLogP_ds = .5 * sum(X.*( -I1 .* d_DD_prob(-s/2,sigma,x) + I2 .* d_DD_prob(+s/2,sigma,x) ) ./ PI, [1,2]);

                % partial derivative wrt d1
                dLogP_dd1 = sum( X.*(dI1 .* DD_prob(-s/2,sigma,x) ./ PI - dI1./(I1+I2)), [1,2]);

                % partial derivative wrt d2
                dLogP_dd2 = sum( X.*(dI2 .* DD_prob(+s/2,sigma,x) ./ PI - dI2./(I1+I2)), [1,2]);

                % gradient of the loglikelihood
                grad = [dLogP_ds, dLogP_dd1, dLogP_dd2]';

                %% LOG-LIKELIHOOD AT CURRENT ESTIMATE 
                loglike(i,j) = sum(X.*log(P),[1,2]);      
    
                %% PERFORM GRADIENT ASCENT UPDATE STEP USING ADAM OPTIMIZER
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
        
        % compute I1, I2
        [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);

        % compute photon arrival prob
        P = (I1 .* DD_prob(-s/2,sigma,x) + I2 .* DD_prob(+s/2,sigma,x)) ./ (I1+I2);

        % store the log likelihood of the current MLE random initialization 
        loglike(i,j+1) = sum(X.*log(P),[1,2]);

    end

    % choose MLE candidate with greatest log likelihood
    [~,k] = max(loglike,[],'all');
    [ki,kj] = ind2sub(size(loglike),k);
    theta_mle = mle(:,ki,kj);

end