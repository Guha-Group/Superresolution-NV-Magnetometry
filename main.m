%%%% SETUP %%%%%

% imager details
sigma = 1;              % Rayliegh limit
integration_time = 100;   % integration time [seconds]
q = 0:10;               % HG mode indices
alpha = .1;             % stage switching parameter for 2-stage receiever alpha \in [0,1];
num_sigma = 10;         % extent of direct imaging detector in number of rayleigh units
num_x = 1001;           % number of detector pixels
x = num_sigma*sigma*linspace(-.5,.5,num_x); % direct imaging detector domain 

% microwave ODMR sweep details
eta0 = 5e4;      % mean photon emission rate of an NV center [photons/s] at zero external magnetic field
chi = .45;       % coupling strength to non-radiative transition (sets NV contrast over ODMR sweep)
linewidth=1;%linewidth = 0.02;% lorentzian width [GHz]
w0 = 0;           %w0 = 2.87;       % zero-field frequency [GHz]
wB = 3*linewidth;   %wB = 3;          % bias-field frequency [GhZ] - 3[Ghz] <-> 4.7 mT bias magnetic field
w_initial=wB-10*linewidth;%w_initial = w0;  % initial freq of ODMR scan [GHz]
w_final=wB+10*linewidth;%w_final = 3.1;   % final freq of ODMR scan [GHz]
omega = (w_initial:linewidth/10:w_final)'; % ODMR microwave freq sweep
schedule = ones(size(omega))/numel(omega);

% setup imaging/spade
IMG.sigma = sigma;
IMG.integration_time = integration_time;
IMG.alpha = alpha;
IMG.q = q;
IMG.x = x;

% setup ODMR
ODMR.eta0 = eta0;
ODMR.chi = chi;
ODMR.w0 = w0;
ODMR.wB = wB;
ODMR.linewidth = linewidth; 
ODMR.schedule = schedule;
ODMR.omega = omega;

%{
% estimation parameters
s = sigma*.25;
d1 = wB - linewidth*.5;
d2 = wB + linewidth*.5;
theta = [s,d1,d2]';

% Test Adaptive Receiver
[theta_mle, pho_tot, alpha_opt, schedule_opt] =...
    Simulate_Adaptive_Receiver(theta,true,IMG,ODMR);

% Test 2-Stage Receiver
[theta_mle,pho_tot] = Simulate_TwoStage_Receiver(theta, IMG, ODMR);

% Test DD-ODMR 
[theta_mle,pho_tot] = Simulate_ODMR_Receiver(theta, IMG, ODMR);
%}

% CRB Limits For 2-stage receiever with adaptive switching
Delta = @(x,beta) beta.*(x/sigma)*linewidth + wB; % Linear ramp detuning
S = sigma*.1;
BETA = 10.^linspace(-1.5,2,1000);

ALPHA_OPT = zeros([numel(S),numel(BETA)]);
CRB_2STAGE = zeros([numel(S),numel(BETA)]);
CRB_ODMR = zeros([numel(S),numel(BETA)]);
for i = 1:numel(S)
    for j = 1:numel(BETA)
        d1 = Delta(-S(i)/2,BETA(j));
        d2 = Delta(+S(i)/2,BETA(j));
        assert(d1<(max(omega)-2*linewidth))
        assert(d1>(min(omega)+2*linewidth))
        assert(d2<(max(omega)-2*linewidth))
        assert(d2>(min(omega)+2*linewidth))
        
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
    plot(S(k).*BETA,ALPHA_OPT(k,:),'LineWidth',1.5)
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
for i = 1:numel(S)
    p = plot(BETA,CRB_2STAGE(i,:),'LineWidth',1.5);
    colors(i,:) = p.Color; 
end
for i = 1:numel(S)
    plot(BETA,CRB_ODMR(i,:),'--','Color',colors(i,:),'LineWidth',1.5);
end
hold off
title('2-Stage Receiver v. DD-ODMR','interpreter','latex')
xlabel({'Detuning Gradient in Units $[w/\sigma]$', '$d\Delta/dx \approx (\Delta_2-\Delta_1)/s$ '},'interpreter','latex')
ylabel({'Cramer-Rao Bound on MSE'},'interpreter','latex')
leg = legend(legend_names,'interpreter','latex');
title(leg,'NV Separation $s$','interpreter','latex');
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([min(BETA),max(BETA)])
axis square
box on


%{
beta = 100;
S = 10.^(linspace(-4,1,1000))*sigma;
ALPHA_OPT = zeros(size(S));
for k = 1:numel(BETA)
    d1 = Delta(-S(k)/2,beta);
    d2 = Delta(+S(k)/2,beta); 
    ALPHA_OPT(k) = OptimizedSwitching(IMG,ODMR,S(k),d1,d2);
end
plot(S,ALPHA_OPT)
set(gca,'XScale','log')
%}
%{
%% parameter scans for Monte-Carlo
s = sigma*10.^linspace(-2,0,10);
d_bar = 3*linewidth; % mean zeeman splitting freq
d_del = linewidth*10.^linspace(-2,0,10); %
d1 = d_bar - d_del;
d2 = d_bar + d_del;
trials = 100;
theta = zeros(numel(s),numel(d_del),3);

for i = 1:numel(s)
    for j = 1:numel(d_del)

        % set the parameter values
        theta(i,j,:) = [s(i),d1(j),d2(j)];

        for t=1:trials
            % simulate estimation for all of the receivers
            receivers(1).theta_mle(i,j,:,t) = Simulate_ODMR_Receiver(theta(i,j,:),IMG,ODMR);
            receivers(2).theta_mle(i,j,:,t) = Simulate_TwoStage_Receiver(theta(i,j,:),IMG,ODMR);
            %receievers(3).theta_mle(i,j,:,t) = Simulate_Adaptive_Receiver(theta(i,j,:),IMG,ODMR);
        end
    end
end
%}
