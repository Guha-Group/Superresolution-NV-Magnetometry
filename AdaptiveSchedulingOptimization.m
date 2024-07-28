function [omega_opt_i,cost_opt] = AdaptiveSchedulingOptimization(IMG,ODMR,theta)
    
    s = theta(1);
    d1 = theta(2);
    d2 = theta(3);

    % compute CFIM of ODMR stage
    cfim = CFIM_ODMR(IMG,ODMR,s,d1,d2);
    
    % transform to gradient coordinates
    J = Jacobian(s,d1,d2);
    cfim = pagemtimes(J.',pagemtimes(cfim,J));
    
    % compute the CRB
    crb = pageinv(cfim);

    % define the cost function
    cost = crb(3,3,:);

    % optimize the cost function with respect to the best detuning frequency
    [cost_opt,omega_opt_i] = max(squeeze(cost));
end