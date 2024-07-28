function dp_dx0 = d_HG_prob(x0,sigma,q)
    % computes the derivative of the Hermite-Gauss modal probability
    % distribution with respect to the position of a point source x0.

    y = (x0/2/sigma).^2;
    
    qq = q-1; qq(qq<0) = 0;
    
    dp_dx0 = (y.^(qq)./factorial(qq).*(q>=1) - y.^(q)./factorial(q)).*exp(-y).*x0/(2*sigma^2);
end