function [alpha_opt,crb] = OptimizedSwitching(IMG,ODMR,s,d1,d2)
    
    % Classical Fisher Info
    CFIM_SP = CFIM_SPADE(IMG,ODMR);
    CFIM_DD = sum(CFIM_ODMR(IMG,ODMR,s,d1,d2),3);
    
    % get jacobian for coordiante transform
    J=Jacobian(s,d1,d2);

    % optimize the switching rate
    alpha = permute(linspace(0,1,1e4),[3,1,2]);

    % total CFIM for params [s,d1,d2]
    CFIM_y = alpha.*CFIM_SP + (1-alpha).*CFIM_DD;
    CFIM_x = pagemtimes(pagetranspose(J),pagemtimes(CFIM_y,J));
    CRB = pageinv(CFIM_x);

    %[~,id] = min(sum(eye(3).*CRB,[1,2]),[],3);
    [crb,id] = min(CRB(3,3,:),[],3);
    alpha_opt = alpha(id);
 
end


