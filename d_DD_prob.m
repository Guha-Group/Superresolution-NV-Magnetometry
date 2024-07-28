function dp_dx0 = d_DD_prob(x0,sigma,x)
    % returns the derivative of the direct imaging probability distribution 
    % with respect to the point source location x0.

    dp_dx0 = 2*DD_prob(x0,sigma,x) .* (x-x0)/(2*sigma^2);
end
