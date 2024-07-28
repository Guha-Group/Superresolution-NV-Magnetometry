function [dI1,dI2] = d_ODMR_2nvFlux(omega, d1, d2, w0, linewidth, eta0, chi)
    % computes the derivative of the flux rate from either NV center with
    % respect to its detuning frequency
    % ------------------------------ 
    % ----------- INPUTS -----------
    % ------------------------------
    % omega             : [W,1] the odmr frequencies
    % d1                : [1,1] Zeeman detuning of NV1 due to local magnetic field
    % d2                : [1,1] Zeeman detuning of NV2 due to local magnetic field
    % w0                : [1,1] Zero B-Field center frequency
    % linewidth         : [1,1] Lorentzian linewidth
    % eta0              : [1,1] mean photon emission rate of each NV in the absence of a magnetic field
    % chi               : [1,1] coupling strength to non-radiative triplet-state transition [0,1] 
    % ------------------------------
    % ----------- OUTPUTS ----------
    % ------------------------------
    % dI1                : [W,1] derivative of emission rate of NV1 at each odmr frequency
    % dI2                : [W,1] derivative of emission rate of NV2 at each odmr frequency

    % Define Lorentzian and its derivative wrt to the detuning
    Lorentzian = @(omega, d) 1 ./ (1 + (omega - w0 - d).^2 / linewidth^2);
    dLorentzian = @(omega, d) 2*(omega - w0 - d) / linewidth^2 .* Lorentzian(omega,d).^2;
    
    % flux rates from either NV source
    dI1 = -eta0 * chi/2 * (-dLorentzian(omega,-d1) + dLorentzian(omega,+d1));
    dI2 = -eta0 * chi/2 * (-dLorentzian(omega,-d2) + dLorentzian(omega,+d2));
end
