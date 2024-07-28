function p = HG_prob(x0,sigma,q)
    % returns modal probability distribution over the Hermite-Gauss modes
    % induced by a point source shifted to location x0 on the image plane.

    x = (x0/2/sigma);
    p = abs((x).^q .* exp(-x.^2/2) ./ sqrt(factorial(q))).^2;
end
