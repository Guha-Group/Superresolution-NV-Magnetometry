function [n1_opt,s_var,s_D_var,s_B_var] = OptimizeFirstSwitching(s,N_pho,n1_samples,sigma)
    % calculate optimal stage 1 and 2 photons given the half-separation
    % 
    % s             : half-separation of the sources
    % N_pho         : total number of photons available for first and second stages
    % n1_samples    : granularity of n1 sampling for search. Minimum switching fraction is constrainted to n1_samples/N_pho
    % sigma         : diffraction limit
    
    assert(n1_samples <= N_pho) % make sure we cannot sub-sample the photons allocated to stage 1
    
    % search range of n1 and n2 allocations
    n1 = round(linspace(1,N_pho,n1_samples));
    n2 = N_pho - n1;

    % calculate ML estimator variance for direct imaging (stage 1)
    s_D_var = (n1.*CFI_DD_s(s,sigma)).^(-1);

    % calculate ML estimator variance for BSPADE (stage 2)
    s_B_var = inf(size(n1));      % container for separation estimator variance under BSPADE measurements
    parfor j=1:(numel(n1)-1)
        s_B_var(j) = BSPADE_Variance(s,n1(j),n2(j),sigma);
    end

    % total variance
    s_var = (1./s_D_var + 1./s_B_var).^(-1);

    % get allocation with minimum estimator variance
    [~,opt_id] = min(s_var);
    n1_opt = n1(opt_id);
end


function s_B_var =  BSPADE_Variance(s,n1,n2,sigma)

    % calculates the variance of the misalignment-averaged ML estimator of
    % the full separation from a BSPADE measurement.
    theta = 2*s;

    % setup the pointing error distribution
    sigma_epsilon = sqrt(sigma^2/n1);
    
    %epsilon = sigma_epsilon*linspace(-4,4,5e2);
    epsilon = sigma_epsilon*linspace(0,4,251);
    de = epsilon(2)-epsilon(1);
    %prob_epsilon = 1/sqrt(2*pi*sigma_epsilon^2).*exp(-(epsilon./sigma_epsilon).^2 / 2);
    prob_epsilon = 2* 1/sqrt(2*pi*sigma_epsilon^2).*exp(-(epsilon./sigma_epsilon).^2 / 2);

    % calculate the probability of photon arrival in the HG00 mode from
    % source 1 and source 2 under different pointing errors
    q1 = ((epsilon-theta/2)/2/sigma).^2;
    q2 = ((epsilon+theta/2)/2/sigma).^2;
    
    % calculate total photon arrival probability in the HG00 mode for each
    % pointint error
    g0 = (1/2)*(exp(-q1)+exp(-q2));

    % calculate the BSPADE PDF under each pointing error
    k = (0:n2)';
    [G0,K] = meshgrid(g0,k);
    N2 = n2*ones(size(K));
    %p_B_Eps = binopdf(K,N2,G0); % conditional distribution of BSPADE given pointing error
    %p_B_Eps = normpdf(K,N2.*G0,sqrt(N2.*G0.*(1-G0))); % approximation of binomial distribution with gaussian

    % partition parameters into those used for binomial pdf
    % and those suitable for the gaussian approximation of a binomial pdf
    binom_sampling_threshold = 10; % binomial sampling threshold for gaussian approximation
    use_binom = (g0 < binom_sampling_threshold/n2) | ((1-g0) < binom_sampling_threshold/n2);
    G0_1 = G0(:,use_binom);  K_1 = K(:,use_binom);  N2_1 = N2(:,use_binom);
    G0_2 = G0(:,~use_binom); K_2 = K(:,~use_binom); N2_2 = N2(:,~use_binom);
    
    % conditional distribution of BSPADE given pointing error
    p_B_Eps = zeros(size(K));
    p_B_Eps(:,use_binom) =  binopdf(K_1,N2_1,G0_1);
    p_B_Eps(:,~use_binom) = normpdf(K_2,N2_2.*G0_2,sqrt(N2_2.*G0_2.*(1-G0_2)));
    
    % marginalize distribution of BSPADE with respect to pointing error 
    p_B = sum(p_B_Eps.*prob_epsilon,2)*de;

    % calculate the conditional MLE for each possible BSPADE measurement outcome
    % given some pointing error
    theta_B = BSPADE_MLE(k,n2,epsilon,sigma);

    % calculate the expected MLE averaged over the pointing error
    % distribution
    expected_theta_B = sum(theta_B .* prob_epsilon, 2)*de;

    % calculate the variance of the expected separation estimator
    theta_B_var = sum( (expected_theta_B - theta).^2 .* p_B);

    % convert to half-separation variance
    s_B_var = theta_B_var/4;
end


function theta_B = BSPADE_MLE(k,n,epsilon,sigma)
    % Michael Grace's separation estimator under centroid misalignment
    % "Approaching quantum-limited imaging resolution without prior
    % knowledge of the object location" -- Grace et. al. 2020

    % checks on input
    assert(all(k<=n))
    assert(n>0)
    assert(all(k>=0))

    % a convenient rate parameter
    Q0 = (epsilon/2/sigma).^2;

    % Solve Cubic Equations in vectorized forms

    % coefficients to cubic Taylor expansion
    a3 = -(1/6) + Q0 - (2/3)*Q0.^2 + (4/45)*Q0.^3;
    a2 = (1/2) - 2*Q0 + (2/3)*Q0.^2;
    a1 = -1 + 2*Q0;
    a0 = 1-(k/n)*exp(Q0);

    % normalized form of cubic equation
    c0 = a0./a3;
    c1 = a1./a3;
    c2 = a2./a3;

    % solve cubic using Wolters Algorithm : https://quarticequations.com/Cubic.pdf
    q = (c1/3) - (c2/3).^2;
    q = repmat(q,[numel(k),1]);
    r = (c1.*c2 - 3*c0)/6 - (c2/3).^3;
    
    % some cases
    one_real = (r.^2 + q.^3)>0;
    rgz = r>=0;
    qgz = q>=0;

    % case 1 solution (only one real solution)
    A = ( abs(r) + (r.^2 + q.^3).^(1/2) ).^(1/3);
    A(~one_real) = 0;
    t1 = A - q./A; t1(~rgz) = -t1(~rgz);
    t1(~one_real) = 0;
    z11 = (t1 - c2/3).*one_real;

    % case 2 solution (three real solutions - pick largest)
    t2 = acos(r./(-q).^(3/2));
    t2(one_real | qgz) = 0;
    phi1 = t2/3;

    
    % All three real solutions with ordering z32 <= z22 <= z12
    %z12 = (2*sqrt(-q).*cos(phi1 + 0) - c2/3).*(~one_real);       % largest
    %z22 = (2*sqrt(-q).*cos(phi1 - 2*pi/3) - c2/3).*(~one_real);  % middle
    z32 = (2*sqrt(-q).*cos(phi1 + 2*pi/3) - c2/3).*(~one_real);   % smallest
    
    % compile the solutions
    mu_B = z11 + z32;   % smallest separation solution

    
    % Run a few iterations of Newton-Raphson to improve precision
    for j = 1:5
        mu_B = mu_B - (mu_B.^3 + c2.*mu_B.^2 + c1.*mu_B + c0)./(3*mu_B.^2 + 2*c2.*mu_B + c1); 
    end   
    
    %{
    % check against numerical solver
    syms x
    i = 300;
    j = 85;
    S = vpasolve(x^3 + c2(j)*x^2 + c1(j)*x + c0(i,j),x)
    one_real(i,j)
    mu_B(i,j)
    numerical_error = mu_B.^3 + c2.*mu_B.^2 + c1.*mu_B + c0;
    %}

    % set negative estimates to zero
    mu_B(mu_B<0) = 0;
     
    % transform 
    theta_B = 4*sigma*sqrt(mu_B);
end

function CFI_1_s = CFI_DD_s(s,sigma)
    % calculates the CFI for the half separation under direct imaging;
    x = sigma*linspace(-10,10,1e4); dx = x(2)-x(1);
    px = exp(-(s/sigma).^2/2)/sqrt(2*pi*sigma^2).*...
              exp(-x.^2/2/sigma^2).*cosh(x.*s/sigma^2);
    dps = ((x.*tanh(x.*s/sigma^2)-s)/sigma^2).^2;
    CFI_1_s = sum(px.*dps).*dx;

end