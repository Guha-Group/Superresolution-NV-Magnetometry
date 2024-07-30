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
    theta = DS.theta(:,array_id); % get parameters for the given array_id
    IMG = DS.IMG;
    ODMR = DS.ODMR;
    adapt_schedule = DS.adapt_schedule;
    
    % make the save directory
    mkdir(DS.save_dir)

    % run parameter scans using matlab's Parallel Computing Toolbox
    %parpool(num_workers)

    for t=1:DS.trials
    %parfor t=1:DS.trials
        
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
                cfg_data(t).alpha_opt = alpha_opt;
                cfg_data(t).schedule_opt = schedule_opt;

        end

        % store the parameter estimates
        cfg_data(t).theta_mle = theta_mle;
        cfg_data(t).pho_tot = pho_tot;
    end

    DS.data(cfg_id{:}) = {cfg_data};
    
    % save current data structure
    fname = [num2str(array_id),'_cfg','.mat'];
    save(fullfile(DS.save_dir,fname),'cfg_id','DS')
    
end