function cfim = CFIM_ODMR(IMG,ODMR,s,d1,d2)
    % calculates the CFIM for the ODMR measurement multiplied by photon
    % rate induced over the microwave drive and the relative schedule time 
    % spent lingering at each microwave drive. The resulting CFIM is in
    % units of information per second.

    % useful constants
    gamma = 1/(2*IMG.sigma)^2;

    % Emitted rates of either NV
    [I1,I2] = ODMR_2nvFlux(ODMR.omega,d1,d2,ODMR.w0,ODMR.linewidth,ODMR.eta0,ODMR.chi);
    I1 = permute(I1,[3,2,1,4,5]);
    I2 = permute(I2,[3,2,1,4,5]);
    eta2 = (I1+I2);

    % compute derivatives of I1 (I2) wrt to detuning d1 (d2)
    [dI1,dI2] = d_ODMR_2nvFlux(ODMR.omega,d1,d2,ODMR.w0,ODMR.linewidth,ODMR.eta0,ODMR.chi);
    dI1 = permute(dI1,[3,2,1,4,5]);
    dI2 = permute(dI2,[3,2,1,4,5]);

    % construct outer product vector for CFIM
    u = cat(1,(I1-I2)./eta2,...
              +2*s.*I2.*dI1./(eta2.^2),...
              -2*s.*I1.*dI2./(eta2.^2));


    cfim = gamma*u.*permute(u,[2,1,3,4,5]);
    cfim = eta2.*permute(ODMR.schedule,[3,2,1]).*cfim;
    %cfim = sum(cfim,3);
end