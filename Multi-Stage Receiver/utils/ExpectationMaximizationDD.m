function [x1,x2,b1,b2] = ExpectationMaximizationDD(X,sigma, max_iterations,position_flag, x12_init, brightness_flag, b12_init)
    % Expectation Maximization algorithm for a gaussian binomial mixture

    % number of photons
    N = numel(X);

    % get some basic initializing estimates
    if position_flag
        x0 = mean(X);
        s = abs(normrnd(0,sigma));
        x1 = x0 - s;
        x2 = x0 + s;
    else
        x1 = x12_init(1);
        x2 = x12_init(2);
    end

    if brightness_flag
        b1 = 0.5;
        b2 = 0.5;
    else
        b1 = b12_init(1);
        b2 = b12_init(2);
    end
    
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
        if position_flag
            x1_next = sum(w1.*X) / sum(w1);
            x2_next = sum(w2.*X) / sum(w2);
        else
            x1_next = x1;
            x2_next = x2;
        end
        
        % update estimates of relative brightness
        if brightness_flag
            b1_next = sum(w1)/N;
            b2_next = sum(w2)/N;
        else
            b1_next = b1;
            b2_next = b2;
        end
        % determine the change in the parameters
        differences = abs([(x1_next - x1),(x2_next - x2),(b1_next - b1),(b2_next - b2)]);

        % update convergence criteria
        k = k+1;
        converged = (k>=max_iterations) || (all(differences < error_epsilon));

        % update estimated positions and brightnesses
        x1 = x1_next;
        x2 = x2_next;
        b1 = b1_next;
        b2 = b2_next;
        
    end

    % flip source labels to be in ascending order with respect to position
    x_est =  [x1,x2]*(x1<=x2) + [x2,x1]*(x2<x1);
    b_est =  [b1,b2]*(x1<=x2) + [b2,b1]*(x2<x1);
    x1 = x_est(1); x2 = x_est(2);
    b1 = b_est(1); b2 = b_est(2);
end