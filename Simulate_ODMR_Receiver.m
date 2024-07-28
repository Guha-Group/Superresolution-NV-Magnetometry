function [theta_mle,pho_tot] = Simulate_ODMR_Receiver(theta, IMG, ODMR)
    % Simulates an instance of the ODMR Receiver
    %
    % theta = [s,d1,d2] : [N,2] the vector of parameters
    % integration_time  : [1,1] total integration time allotted
    % omega             : [W,1] the odmr frequencies
    % schedule          : [W,1] the odmr schedule (must sum to 1)
    % x                 : [1,M] the direct imaging detector positions
    % eta0              : [1,1] mean photon emission rate of each NV in the absence of a magnetic field 

    % unpack parameters
    s  = theta(1);
    d1 = theta(2);
    d2 = theta(3);
    
    % setup imaging
    x = IMG.x;
    sigma = IMG.sigma;
    integration_time = IMG.integration_time;

    % setup ODMR
    eta0 = ODMR.eta0;
    chi = ODMR.chi;
    w0 = ODMR.w0;
    linewidth = ODMR.linewidth;
    schedule = ODMR.schedule;
    omega = ODMR.omega;

    
    % Get mean photon rates from each NV center as a function of the ODMR
    % scan frequency    
    [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);

    % Get time spent at each frequency through the ODMR schedule
    t = integration_time*schedule;

    % Define the relative brightness of either NV center
    eta = I1 + I2;
    b1 = I1./eta;
    b2 = I2./eta;

    % Get photon probability distribution 
    p1 = DD_prob(-s/2,sigma,x); p1 = p1/sum(p1);
    p2 = DD_prob(+s/2,sigma,x); p2 = p2/sum(p2);
    P_DD = b1.*p1 + b2.*p2;
    
    % Simulate and ODMR measurement (Poisson photon statistics)
    X = SimulateMeasurement(P_DD,eta.*t);

    % get maximimum likelihood estimator
    theta_mle = ODMR_MLE(X, x, IMG, ODMR);

    pho_tot = sum(X,'all');
end

