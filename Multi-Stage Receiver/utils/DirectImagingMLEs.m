function [x0_est, s_est] = DirectImagingMLEs(X,sigma)
    % Calculates the maximum likelihod estimators for the midpoint and the
    % separation of two balanced point sources from a direct imaging
    % measurement captured with a system characterized by a gaussian PSF.
    
    % ML estimator equations
    x0_est = mean(X);
    s_ML_transcendental = @(s_ML, x0) s_ML - mean( tanh((X-x0).*s_ML./sigma^2).*(X-x0),2);
    x0_ML_transcendental = @(x0_ML, s) x0_ML - mean(X + s * tanh((X-x0_ML) .* s / (sigma^2)),2);
    
    % alternate estimators twice to get convergence
    for k = 1:2
        Ls = @(s_ML) s_ML_transcendental(s_ML, x0_est);
        s_est = abs(fzero(Ls,sigma));
        Lx = @(x0_ML) x0_ML_transcendental(x0_ML,s_est);
        x0_est = fzero(Lx,0);
    end
end
