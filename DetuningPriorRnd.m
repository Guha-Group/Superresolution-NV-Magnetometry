function [d1,d2]=DetuningPriorRnd(ODMR,num_samples)
    % Random sampler of the detunings from a user-defined prior. 
    d_diff = abs(normrnd(0,ODMR.linewidth,[1,1,1,1,num_samples]));
    d1 = (ODMR.wB - ODMR.w0) - d_diff/2;
    d2 = (ODMR.wB - ODMR.w0) + d_diff/2;
end