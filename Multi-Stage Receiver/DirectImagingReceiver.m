function params_out = DirectImagingReceiver(x0,s,kappa,N,sigma,max_iterations)
    addpath('utils/')
 
    % Scene Parameters
    % x0 - geometric midpoint of the point sources
    % s  - half-separation of the point sources (must be greater than 0)
    % kappa - relative brightness bias (must be between [-0.5,+0.5]

    % System Parameters
    % sigma - diffraction limit
    % N     - photon budget

    % Absolute Positions
    x1 = x0-s;
    x2 = x0+s;
    
    % Absolute Brightnesses
    b1 = 0.5 - kappa;
    b2 = 0.5 + kappa;
   
    % Simulate a direct imaging measurement
    X = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],N,sigma);
    
    % Estimate source brightnesses via expectation maximization
    brightness_flag = 1;
    position_flag = 1;
    [x1_est, x2_est, b1_est, b2_est] = ExpectationMaximizationDD(X,sigma,max_iterations,position_flag,[nan,nan],brightness_flag,[nan,nan]);

    % map to target parameter estimates
    x0_est = (x1_est + x2_est)/2;
    s_est = (x2_est - x1_est)/2;
    kappa_est = (b2_est-b1_est)/2;

    % collect parameter estimates
    params_out = [x0_est; s_est; kappa_est];
end
