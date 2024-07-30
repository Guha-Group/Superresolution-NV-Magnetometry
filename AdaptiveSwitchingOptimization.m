function [alpha_opt,cost_opt] = AdaptiveSwitchingOptimization(IMG, ODMR, s)
    % Computes the optimal switching time given an intermediate estimate of
    % the separation parameter
    % IMG       : imaging system struct
    % ODMR      : ODMR system struct
    % s         : NV separation
    
    % get first stage CFIM
    cfim_spade = CFIM_SPADE(IMG,ODMR);

    % get expected value of cost function
    num_samples = 1e4;
    [d1_samples,d2_samples] = DetuningPriorRnd(ODMR, num_samples); 
    s = s*ones(size(d1_samples)); % extend dimensions of s
    
    % cfim of odmr for all samples of the prior
    cfim_odmr = CFIM_ODMR(IMG,ODMR,s,d1_samples,d2_samples);
    cfim_odmr = sum(cfim_odmr,3);

    % samples of switching parameter
    alpha = permute(linspace(1e-6,1-1e-6,1e3),[1,3,4,2]);

    % get two-stage cfim
    cfim_2stage = alpha.*cfim_spade + (1-alpha).*cfim_odmr;

    % apply Jacobian transformation
    J = Jacobian(s,d1_samples,d2_samples);
    cfim_2stage = pagemtimes(pagetranspose(J),pagemtimes(cfim_2stage,J));
    
    % get crbs for all samples
    crb_2stage = pageinv(cfim_2stage);
    
    % expected crb over delta prior
    crb_expectation = mean(crb_2stage,5);

    % return the switching time with the highest expected information
    %cost = sum(eye(3).*crb_expectation,[1,2]);
    cost = crb_expectation(3,3,1,:);
    [cost_opt,id] = min(cost,[],4);
    alpha_opt = alpha(id); 
end