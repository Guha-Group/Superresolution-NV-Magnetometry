function photon_data = SimulateMeasurement(P,U)
% P: [W,M] array of probability distributions. 
%           Each row is a discrete probability distribution.
% U: [W,1] array of mean numbers of photons collected per ODRM frequency.
%
%    Dim 1 (W) is the ODMR frequency axis
%    Dim 2 (M) is the mode index axis
    photon_data = poissrnd(sum(P.*U,3));
end
