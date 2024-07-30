%%%% SETUP %%%%%

% imager details
sigma = 1;              % Rayliegh limit
integration_time = 1;   % integration time [seconds]
q = 0:10;               % HG mode indices
alpha = .1;             % stage switching parameter for 2-stage receiever alpha \in [0,1];
num_sigma = 10;         % extent of direct imaging detector in number of rayleigh units
num_x = 1001;           % number of detector pixels
x = num_sigma*sigma*linspace(-.5,.5,num_x); % direct imaging detector domain 

% microwave ODMR sweep details
eta0 = 5e4;             % mean photon emission rate of an NV center [photons/s] at zero external magnetic field
chi = .45;              % coupling strength to non-radiative transition (sets NV contrast over ODMR sweep)
linewidth = 0.02;       % lorentzian width [GHz]
w0 = 2.87;              % zero-field frequency [GHz]
wB = w0+5*linewidth;    % bias-field frequency [GhZ] - 3[Ghz] <-> 4.7 mT bias magnetic field
dB = wB-w0;             % detuning bias displacement
w_initial = max(w0, wB - 5*linewidth);      % initial freq of ODMR scan [GHz]
w_final = wB + 5*linewidth;                 % final freq of ODMR scan [GHz]     
omega = (w_initial:linewidth/10:w_final)';  % ODMR microwave freq sweep
schedule = ones(size(omega))/numel(omega);  % ODMR microwave drive schedule

% ADAM optimizer gradient ascent details for MLE
num_inits = 30;     % number of times to initalize the gradient ascent algorithm (robustness to local maxima)
num_steps = 500;    % number of gradient ascent iterations
learning_rate = 1.5/eta0;  % step size in the gradient.
beta_1 =  .9;
beta_2 = .99;
epsilon = 1e-8;

% setup imaging/spade
IMG.sigma = sigma;
IMG.integration_time = integration_time;
IMG.alpha = alpha;
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


%% Test one-off run of receivers
s  = 0.1*sigma;
d1 = dB;
d2 = dB + 2*linewidth;
theta = [s,d1,d2]';
[theta_mle_1,~,alpha_opt,schedule_opt] = Simulate_Adaptive_Receiver(theta,0,IMG,ODMR);
theta_mle_2 = Simulate_ODMR_Receiver(theta,IMG,ODMR);

%% Optimal Switching analysis for 2-stage gradient sensing

% CRB Limits For 2-stage receiever with adaptive switching
Delta = @(x,beta) beta .* (x/sigma) * linewidth + dB; % Linear ramp detuning
S = sigma*linspace(.01,.5,5)';
beta_limit = ((dB-linewidth)/linewidth * 2*sigma./S) ;
%beta_limit = min(abs(w_initial-delta_B),abs(w_final-delta_B)) / linewidth * 2 * sigma ./ S;
BETA = beta_limit*linspace(1e-3,1,50);
%BETA= beta_limit*10.^linspace(-5,0,50); % magnetic field gradient scaling

ALPHA_OPT = zeros(size(BETA));
CRB_2STAGE = zeros(size(BETA));
CRB_ODMR = zeros(size(BETA));
for i = 1:numel(S)
    for j = 1:size(BETA,2)
        d1 = Delta(-S(i)/2,BETA(i,j));
        d2 = Delta(+S(i)/2,BETA(i,j));
        assert((d1+w0)<=(max(omega)))
        assert((d1+w0)>=(min(omega)))
        assert((d2+w0)<=(max(omega)))
        assert((d2+w0)>=(min(omega)))
        
        % Get optimal 2-Stage switching parameter and corresponding crb for
        % magnetic field gradient estimation
        [a_opt,c_2stage] = OptimizedSwitching(IMG,ODMR,S(i),d1,d2);
        
        % get standard ODMR CRB for magnetic field gradient estimation
        cfim_odmr = sum(CFIM_ODMR(IMG,ODMR,S(i),d1,d2),3);
        J = Jacobian(S(i),d1,d2);
        cfim_odmr = J.'*cfim_odmr*J;
        crb_odmr = inv(cfim_odmr);
        c_odmr = crb_odmr(3,3);

        % collect values
        ALPHA_OPT(i,j) = a_opt;
        CRB_2STAGE(i,j) = c_2stage;
        CRB_ODMR(i,j) = c_odmr;
    end
end


%% PLOTS
figure
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
nexttile(1)
hold on
for k = 1:numel(S)
    %plot(S(k)*BETA(k,:),ALPHA_OPT(k,:),'LineWidth',1.5) % alpha opt is a function of the detuning difference only
    plot(BETA(k,:),ALPHA_OPT(k,:),'LineWidth',1.5)    
end
hold off
xlabel({'Detuning Gradient in Units $[w/\sigma]$', '$d\Delta/dx \approx (\Delta_2-\Delta_1)/s$ '},'interpreter','latex')
ylabel({'Optimal Fraction of Time Allocated to SPADE','$\alpha^*$'},'interpreter','latex')
title('Optimal Stage Switching','interpreter','latex')
%set(gca,'yscale','log')
%set(gca,'xscale','log')
legend_names = arrayfun(@(k) sprintf('$%.2f\\sigma$',S(k)/sigma),1:numel(S),'UniformOutput',false);
leg = legend(legend_names,'interpreter','latex');
title(leg,'NV Separation $s$','interpreter','latex');
axis square
box on

nexttile(2)
colors = zeros(numel(S),3);
hold on
for k = 1:numel(S)
    %p = plot(S(k)*BETA(k,:),S(k)^4*CRB_2STAGE(k,:),'LineWidth',1.5); %%CRBs are globally modulated by factor of s^4
    p = plot(S(k)*BETA(k,:),CRB_2STAGE(k,:),'LineWidth',1.5);
    colors(k,:) = p.Color; 
end
for k = 1:numel(S)
    %plot(S(k)*BETA(k,:),S(k)^4*CRB_ODMR(k,:),'--','Color',colors(k,:),'LineWidth',1.5);
    plot(S(k)*BETA(k,:),CRB_ODMR(k,:),'--','Color',colors(k,:),'LineWidth',1.5);
end
hold off
title('2-Stage Receiver v. DD-ODMR','interpreter','latex')
xlabel({'Detuning Gradient in Units $[w/\sigma]$', '$d\Delta/dx \approx (\Delta_2-\Delta_1)/s$ '},'interpreter','latex')
ylabel({'Cramer-Rao Bound on MSE'},'interpreter','latex')
leg = legend(legend_names,'interpreter','latex');
title(leg,'NV Separation $s$','interpreter','latex');
set(gca,'yscale','log')
set(gca,'xscale','log')
axis square
box on


