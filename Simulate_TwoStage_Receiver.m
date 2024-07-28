function [theta_mle,pho_tot] = Simulate_TwoStage_Receiver(theta, IMG, ODMR)
% Here we simulate a 2-stage receiever for 

    % unpack parameters
    s  = theta(1);
    d1 = theta(2);
    d2 = theta(3);
    
    % setup imaging
    sigma = IMG.sigma;
    integration_time = IMG.integration_time;
    alpha = IMG.alpha;
    x = IMG.x;
    q = IMG.q;

    % setup ODMR
    eta0 = ODMR.eta0;
    chi = ODMR.chi;
    w0 = ODMR.w0;
    linewidth = ODMR.linewidth;
    schedule = ODMR.schedule;
    omega = ODMR.omega;

    % Get timing schedules for each stage
    t1 = alpha * integration_time;
    t2 = (1-alpha) * integration_time * schedule;
    
    %% STAGE 1: SPADE
    % Get HG mode probability distribution
    P1 = 0.5*(HG_prob(-s/2,sigma,q) + HG_prob(+s/2,sigma,q));
    Q = SimulateMeasurement(P1, 2*eta0*t1);


    %% STAGE 2: ODMR
    % Get mean photon rates from each NV center as a function of the ODMR
    % scan frequency    
    [I1,I2] = ODMR_2nvFlux(omega,d1,d2,w0,linewidth,eta0,chi);

    % Define the relative brightness of either NV center
    eta = I1 + I2;
    b1 = I1./eta;
    b2 = I2./eta;

    % Get direct imaging ODMR probability distribution 
    p1 = DD_prob(-s/2,sigma,x); p1 = p1/sum(p1);
    p2 = DD_prob(+s/2,sigma,x); p2 = p2/sum(p2);
    P2 = b1.*p1 + b2.*p2;

    X = SimulateMeasurement(P2, eta.*t2);


    %% ESTIMATION
    theta_mle = TWOSTAGE_MLE(Q, q, X, x, IMG,ODMR);

    pho_tot = sum(Q,'all') + sum(X,'all');

end
