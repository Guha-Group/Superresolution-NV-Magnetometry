% Brightness Estimation from Direct Imaging

% System Parameters
sigma = 1;          % diffraction limit
x0 = 0;             % geometric midpoint of the point sources
s = 0.3*sigma;      % half-separation of the point sources (must be greater than 0)
kappa = 0.3;        % relative brightness bias (must be between [-0.5,+0.5]
mu = (s/2/sigma)^2;

% Absolute Positions
x1 = x0-s;
x2 = x0+s;

% Absolute Brightnesses
b1 = 0.5 - kappa;
b2 = 0.5 + kappa;

% Measurement Parameters
N = 1e5;        % total photon budget

% Simulate a direct imaging measurement
X = SimulateDirectImagingMeasurement([x1,x2],[b1,b2],N,sigma);

% Estimate source brightnesses via expectation maximization
max_iterations = 1000;
[x1_est,x2_est,b1_est,b2_est] = ExpectationMaximizationDD(X,sigma,max_iterations);


function [x1,x2,b1,b2] = ExpectationMaximizationDD(X,sigma, max_iterations)
    % Expectation Maximization algorithm for a gaussian binomial mixture

    % number of photons
    N = numel(X);

    % get some basic initializing estimates
    x0 = mean(X);
    s = abs(normrnd(0,sigma));

    % initialize parameters
    x1 = x0 - s;
    x2 = x0 + s;
    b1 = 0.5;
    b2 = 0.5;
    
    % a small number
    error_epsilon = 1e-4;

    % run expectation maximization for finite number of iterations
    converged = 0;
    k = 0;
    while ~converged
        
        % sample weights
        ps1 = normpdf(X,x1,sigma);
        ps2 = normpdf(X,x2,sigma);
        w1 = b1*ps1./(b1*ps1 + b2*ps2);
        w2 = b2*ps2./(b1*ps1 + b2*ps2);

        % update estimates of positions 
        x1_next = sum(w1.*X) / sum(w1);
        x2_next = sum(w2.*X) / sum(w2);
        
        % update estimates of relative brightness
        b1_next = sum(w1)/N;
        b2_next = sum(w2)/N;

        % determine the change in the parameters
        differences = abs([(x1_next - x1),(x2_next - x2),(b1_next - b1),(b2_next - b2)]);

        % update convergence criteria
        k = k+1;
        converged = (k>=max_iterations) || (all(differences < error_epsilon));
        
    end
end