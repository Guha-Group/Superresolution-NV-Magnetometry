function p = DD_prob(x0,sigma,x)
    % returns the direct imaging probability distrubition induced by a
    % point source shifted to loaction x0 on the image plane.
    p = (2*pi*sigma^2).^(-1/2) .* exp(- (x-x0).^2/(2*sigma^2));
end
