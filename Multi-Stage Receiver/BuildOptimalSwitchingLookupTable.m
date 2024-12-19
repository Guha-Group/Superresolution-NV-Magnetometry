% Parameter setups
sigma = 1;
n1_samples = 1000;
s_range = sigma*linspace(1e-2,1,100);
N_range = [1e3,1e4,1e5]';
N1_opt = zeros(numel(N_range), numel(s_range));
S_Var = zeros([size(N1_opt), n1_samples]);
S_D_Var = S_Var;
S_B_Var = S_Var;

for n = 1:numel(N_range)
    N = N_range(n);
    for k = 1:numel(s_range)
        s = s_range(k);
        
        % calculate optimal switching times
        [n1_opt,s_var,s_D_var,s_B_var] = OptimizeFirstSwitching(s,N,n1_samples,sigma);
        
        % assign results
        N1_opt(n,k) = n1_opt;
        S_Var(n,k,:) = s_var;
        S_D_Var(n,k,:) = s_D_var;
        S_B_Var(n,k,:) = s_B_var;
    end
end
opt_switch_frac = N1_opt./N_range;

save('SwitchingLookupTables.mat','opt_switch_frac','N_range','s_range','S_Var','S_D_Var','S_B_Var')
