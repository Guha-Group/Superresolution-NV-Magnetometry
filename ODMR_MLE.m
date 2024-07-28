function theta_mle = ODMR_MLE(X, x, IMG, ODMR)
    % compute MLE from the simulated ODMR measurement via gradient ascent
    % on the likelihood
    
    % unpack vars
    sigma = IMG.sigma;
    eta0 = ODMR.eta0;
    chi = ODMR.eta0;
    omega = ODMR.omega;
    w0 = ODMR.w0;
    linewidth = ODMR.linewidth;
    
    num_inits = 30;     % number of times to initalize the gradient ascent algorithm (robustness to local maxima)
    num_steps = 100;    % number of gradient ascent iterations
    learning_rate = 1/eta0;  % step size in the gradient.

    mle = zeros([3,num_inits]); 
    loglike = zeros(num_inits,num_steps+1);
    for i = 1:num_inits
        theta_mle = [sigma*abs(normrnd(0,1,1)), w0 + linewidth*abs(normrnd(0,5,1,2))]';
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
                
                % get log likelihood
                loglike(i,j) = sum(X.*log(P),[1,2]);
                
                % partial derivative wrt s
                dLogP_ds = .5 * sum(X.*( I1 .* d_DD_prob(-s/2,sigma,x) + I2 .* d_DD_prob(+s/2,sigma,x) ) ./ PI, [1,2]);

                % partial derivative wrt d1
                dLogP_dd1 = sum( X.*(dI1 .* DD_prob(-s/2,sigma,x) ./ PI - dI1./(I1+I2)), [1,2]);

                % partial derivative wrt d2
                dLogP_dd2 = sum( X.*(dI2 .* DD_prob(+s/2,sigma,x) ./ PI - dI2./(I1+I2)), [1,2]);

                % gradient
                grad = [dLogP_ds, dLogP_dd1, dLogP_dd2]';

                theta_mle = theta_mle + learning_rate*grad;      
        end
        
        % store the MLE for the current random initializaion
        mle(:,i) = theta_mle;

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
    [~,k] = max(loglike(:,end));
    theta_mle = mle(:,k);
end