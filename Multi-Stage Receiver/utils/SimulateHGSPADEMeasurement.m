function Q = SimulateHGSPADEMeasurement(xs,bs,num_pho,sigma)
    num_sources = numel(bs);    % total number of sources
    ns = round(bs*num_pho);     % number of photons from each source

    Q = zeros(1,num_pho);       % container for photon mode indices
    kk=1;                       % index counter
    for k = 1:num_sources
        lambda = (xs(k)/2/sigma)^2;
        Q(kk:(kk-1+ns(k))) = poissrnd(lambda,1,ns(k));
        kk = kk + ns(k);
    end
end