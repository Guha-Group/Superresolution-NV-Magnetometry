function cfim = CFIM_SPADE(IMG,ODMR)
    % Calculates the CFIM for the SPADE measurement multplied by the constant
    % photon emission rate of the NVs. The resulting CFIM is in units of
    % informationa per second.

    % useful constants
    gamma = 1/(2*IMG.sigma)^2;

    eta1 = 2*ODMR.eta0; % photon arrival rate in first stage
    cfim = zeros(3);
    cfim(1,1) = gamma;
    cfim = eta1*cfim;
    
end
