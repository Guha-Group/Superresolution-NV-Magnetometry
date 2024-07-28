function [d1,d2]=DetuningPriorRnd(ODMR,num_samples)
    
    d_diff = abs(normrnd(0,ODMR.linewidth,[1,1,1,1,num_samples]));
    d1 = ODMR.wB - d_diff/2;
    d2 = ODMR.wB + d_diff/2;

end