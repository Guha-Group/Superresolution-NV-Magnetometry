function DS = DSformat()
    % consruct a data object format for running and storing Monte-Carlos
    
    % save directory
    save_dir = fullfile('Data','ReceiverMonteCarlo_BGradient');

    % imager details
    sigma = 1;              % Rayliegh limit
    integration_time = 10;  % integration time [seconds]
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
    wB = 3;          % bias-field frequency [GhZ] - 3[Ghz] <-> 4.7 mT bias magnetic field IF NOT BIASED, SET wB = w0
    delta_B = wB-w0; % bias offset from zero-field splitting
    w_initial = max(w0, wB - 5*linewidth); % initial freq of ODMR scan [GHz]
    w_final = wB + 5*linewidth;            % final freq of ODMR scan [GHz]  
    omega = (w_initial:linewidth/10:w_final)'; % ODMR microwave freq sweep
    schedule = ones(size(omega))/numel(omega); % timing for each freq.

    % ADAM optimizer gradient ascent details for MLE computation
    num_inits = 30;             % number of times to initalize the gradient ascent algorithm (robustness to local maxima)
    num_steps = 500;            % number of gradient ascent iterations
    learning_rate = 1.5/eta0;   % step size in the gradient.
    beta_1 =  .9;               % gradient momentum parameter 1
    beta_2 = .99;               % gradient momentum parameter 2
    epsilon = 1e-8;             % a small constant

    %%%%%%%% MONTE-CARLO PARAMETER SWEEPS %%%%%%%%
    receiver_names={'ODMR','Static 2-Stage','Adaptive 2-Stage'};
    trials = 1000;     % number of Monte-carlo trials per configuration
    s = .1*sigma;      % NV separation
    beta_limit = min(abs(w_initial-wB)-linewidth,abs(w_final-wB)-linewidth)/linewidth * 2*sigma/s;
    beta = beta_limit*10.^linspace(-2,0,50); % magnetic field gradient scaling
    Delta = @(x,beta) beta.*(x/sigma)*linewidth + delta_B; % Linear ramp detuning
    d1 = Delta(-s/2,beta);
    d2 = Delta(+s/2,beta);
    theta = cat(1,s*ones(size(beta)),d1,d2);


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