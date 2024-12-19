function DS = DSformat()
    % consruct a data object format for running and storing Monte-Carlos
    
    % save directory
    save_dir = fullfile('Data','BGradientMC_Aug8');

    % imager details
    sigma = 1;              % Rayliegh limit
    integration_time = 50;  % integration time [seconds]
    q = 0:10;               % HG mode indices
    alpha = .1;             % stage switching parameter for 2-stage receiever alpha \in [0,1];
    num_sigma = 10;         % extent of direct imaging detector in number of rayleigh units
    num_x = 1001;           % number of detector pixels
    x = num_sigma*sigma*linspace(-.5,.5,num_x); % direct imaging detector domain 

    % microwave ODMR sweep details
    eta0 = 5e4;      % mean photon emission rate of an NV center [photons/s] at zero external magnetic field
    chi = .45;       % coupling strength to non-radiative transition (sets NV contrast over ODMR sweep)
    linewidth = 0.02;% lorentzian width [GHz]
    w0 = 2.87;       % zero-field frequency [GHz]

    % ADAM optimizer gradient ascent details for MLE computation
    num_inits = 30;             % number of times to initalize the gradient ascent algorithm (robustness to local maxima)
    num_steps = 500;            % number of gradient ascent iterations
    learning_rate = 1.5/eta0;   % step size in the gradient.
    beta_1 =  .9;               % gradient momentum parameter 1
    beta_2 = .99;               % gradient momentum parameter 2
    epsilon = 1e-8;             % a small constant

    %%%%%%%% MONTE-CARLO PARAMETER SWEEPS %%%%%%%%
    receiver_names={'ODMR','Static 2-Stage','Adaptive 2-Stage','Adaptive Scheduling'};
    trials = 1000;     % number of Monte-carlo trials per configuration
    s = .01*sigma;     % NV separation
    beta = 10.^linspace(-2,2,50);                   % magnetic field gradient scaling
    delta_B = linewidth*(7 + max(beta).*s/sigma/2); % bias offset from zero-field splitting [GHz]
    wB = w0 + delta_B;                              % absolute bias frequency [GHz]
    Delta = @(x,beta) beta.*x*(linewidth/sigma) + delta_B;  % Linear ramp detuning
    d1 = Delta(-s/2,beta);                          % spin resonance frequencies of NV1
    d2 = Delta(+s/2,beta);                          % spin resonance frequencies of NV2
    theta = cat(1,s*ones(size(beta)),d1,d2);        % all parameters
    
    % setup ODMR drive frequency range
    w_initial = w0 + min(d1)-linewidth;                 % initial freq of ODMR scan [GHz]
    w_final = w0 + max(d2)+linewidth;                     % final freq of ODMR scan [GHz]  
    omega = (w_initial:linewidth*min(beta):w_final)';  % ODMR microwave freq sweep [GHz]
    schedule = ones(size(omega))/numel(omega);  % scheduling time spent at each ODMR freq. [fraction of total integration time]
    
    %%%%%%%%%%%    SETUP STRUCTS     %%%%%%%%%%%%
    % setup imaging/spade
    IMG.sigma = sigma;
    IMG.integration_time = integration_time;
    IMG.alpha = alpha;  % for Two-Stage Receiver Only
    IMG.q = q;
    IMG.x = x;
    IMG.num_inits = num_inits;
    IMG.num_steps = num_steps;
    IMG.learning_rate = learning_rate;
    IMG.beta_1 = beta_1;
    IMG.beta_2 = beta_2;
    IMG.epsilon = epsilon;

    % setup ODMR
    ODMR.eta0 = eta0;
    ODMR.chi = chi;
    ODMR.w0 = w0;
    ODMR.wB = wB;
    ODMR.linewidth = linewidth; 
    ODMR.schedule = schedule;
    ODMR.omega = omega;

    % data structure with some limited functionality
    DS = struct();
    DS.IMG = IMG;
    DS.ODMR = ODMR;
    DS.save_dir = save_dir;
    DS.trials = trials;
    DS.theta = theta;
    DS.beta = beta;
    DS.adapt_schedule = false;
    DS.receiver_names = receiver_names;
    DS.cfg_size = [numel(receiver_names),numel(beta),numel(s)]; % the dimensionality of the parameter space range
    DS.data = cell(DS.cfg_size);

end