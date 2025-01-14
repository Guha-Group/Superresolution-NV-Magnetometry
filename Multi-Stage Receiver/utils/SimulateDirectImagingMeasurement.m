function X = SimulateDirectImagingMeasurement(xs,bs,num_pho,sigma)
    % Description: Returns the arrival location of photons on an ideal
    % photo-detector.
    num_sources = numel(bs);    % total number of sources
    ns = round(bs*num_pho);     % number of photons from each source
    
    X = zeros(1,num_pho);         % container for photon arrival locations
    kk=1;                       % index counter
    for k = 1:num_sources
        X(kk:(kk-1+ns(k))) = sigma*randn(1,ns(k))+xs(k);
        kk = kk+ns(k);
    end
end