function ReceiverSurvey(array_id, num_workers)
    % runs a monte-carlo survey of ODMR Receiver performance 
    % (adapted for deployment on a compute cluster)
    
    % get the job array id
    if ischar(array_id)
        array_id = str2double(array_id);
    end
    
    % make the DS structure
    DS = DSformat();

    % get configuration indices
    [r,b,s] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {r,b,s};

    % unpack parameters and flags
    theta = DS.theta(:,b,s); % get parameters for the given array_id
    IMG = DS.IMG;
    ODMR = DS.ODMR;
    adapt_schedule = DS.adapt_schedule;
    
    % make the save directory
    mkdir(DS.save_dir)
    
    % containers for simulation outputs
    cfg_theta_mle = zeros(numel(theta),DS.trials);
    cfg_pho_tot = zeros(1,DS.trials);
    cfg_alpha_opt = zeros(1,DS.trials);
    cfg_schedule_opt = zeros(numel(DS.ODMR.omega),DS.trials);

    % run parameter scans using matlab's Parallel Computing Toolbox
    parpool(num_workers)

    %for t=1:DS.trials
    parfor t=1:DS.trials
        
        % set the random number generator seed (must include the parallel
        % worker in the seed for different instances to develop different scenes)
        rng(array_id + t)
        
        switch r
            case 1
                [theta_mle,pho_tot] = Simulate_ODMR_Receiver(theta,IMG,ODMR);
            case 2
                [theta_mle,pho_tot] = Simulate_TwoStage_Receiver(theta,IMG,ODMR);
            case 3
                [theta_mle,pho_tot,...
                alpha_opt, schedule_opt] = Simulate_Adaptive_Receiver(theta, adapt_schedule, IMG,ODMR);
                cfg_alpha_opt(1,t) = alpha_opt;
                cfg_schedule_opt(:,t) = schedule_opt;

        end

        % store the parameter estimates
        cfg_theta_mle(:,t) = theta_mle;
        cfg_pho_tot(1,t) = pho_tot;
    end

    cfg_data.theta_mle = cfg_theta_mle;
    cfg_data.pho_tot = cfg_pho_tot;
    cfg_data.alpha_opt = cfg_alpha_opt;
    cfg_data.schedule_opt = cfg_schedule_opt;

    DS.data(cfg_id{:}) = {cfg_data};
    
    % save current data structure
    fname = [num2str(array_id),'_cfg','.mat'];
    save(fullfile(DS.save_dir,fname),'cfg_id','DS')
    
end