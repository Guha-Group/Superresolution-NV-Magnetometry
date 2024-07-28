function theta_mle = TWOSTAGE_MLE(Q, q, X, x,IMG,ODMR)
    % Q     : [1,N] measurement vector containing the number of photon arrivals in each mode of the SPADE measurement
    % q     : [1,N] vector of mode indices
    % X     : [W,M] measurement matrix containing number of photon arrival at each focal plane position for each odmr frequency
    % x     : [1,M] the direct imaging detector positions

    sigma = IMG.sigma;
    omega = ODMR.omega;
    w0 = ODMR.w0;
    linewidth = ODMR.linewidth;
    eta0 = ODMR.eta0;
    chi = ODMR.chi;

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
            
            %% STAGE 1 GRADIENT:

            % probability over sample space
            P1 = 0.5*(HG_prob(-s/2,sigma,q) + HG_prob(+s/2,sigma,q));
            
            % partial derivative wrt to s
            dLogP1_ds = .5* .5 * sum(Q .* (-d_HG_prob(-s/2,sigma,q) + d_HG_prob(+s/2,sigma,q)) ./ P1 ,2);

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
            dLogP_ds = .5 * sum(X.*( I1 .* d_DD_prob(-s/2,sigma,x) + I2 .* d_DD_prob(+s/2,sigma,x) ) ./ PI, [1,2]);

            % partial derivative wrt d1
            dLogP_dd1 = sum( X.*(dI1 .* DD_prob(-s/2,sigma,x) ./ PI - dI1./(I1+I2)), [1,2]);

            % partial derivative wrt d2
            dLogP_dd2 = sum( X.*(dI2 .* DD_prob(+s/2,sigma,x) ./ PI - dI2./(I1+I2)), [1,2]);

            % gradient
            grad2 = [dLogP_ds, dLogP_dd1, dLogP_dd2]';
            

            %% LOG-LIKELIHOOD AT CURRENT ESTIMATE 
            loglike(i,j) = sum(Q.*log(P1),2) + sum(X.*log(P2),[1,2]);

            %% PERFORM GRADIENT ASCENT UPDATE STEP
            grad = grad1 + grad2;
            theta_mle = theta_mle + learning_rate*grad;


        end
        % store the MLE for the current random initializaion
        mle(:,i) = theta_mle;

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
    [~,k] = max(loglike(:,end));
    theta_mle = mle(:,k);

end

