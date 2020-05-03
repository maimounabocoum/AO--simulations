function [eta,Ppth,Grt,s] = LaserPerformance(s)
% Function to perform analytical steady state analysis of linear 
% or ring lasers calculating:
%   - Slope efficiency
%   - Lasing threshold
%   - Lasing wavelength
%   - Round trip gain of laser
%   - Optimal fiber length
% 
% Inputs:
%
% Outputs: 
%   eta - efficiency of conversion of pump energy to lasing energy, above
%       threshold. given in decimal ratio, not percentage
%   Ppth - pump power at lasing threshold. given in watts.
%   Grt - round trip gain required to lase. given as linear value, not dB.
%  
% Comments:
%   - Uses some quasi-two level assumptions for ytterbium, as noted in 
%       definitions section

% References: 
%   Pfieffer - 1533 nm EDFL as in "Power Characteristic of Erbium-Doped
%       Lasers", IEEE Ph Tech Letters 1992, vol 4 no 8 pp 847-849 (stacks)
%   Barnard - 1533 nm EDFL as in "Rare Earth Doped Fiber Amplifiers and
%       Lasers", IEEE J. Q. Electronics 1994,vol30,no8,pp1817-1830 (stks)
%   Poulsen - 1560 nm EDFL as in "Highly Optimized Er3+ Doped Fiber Ring
%       Laser", IEEE Ph Tech Letters 1993, vol 5 no 6, pp 646-648

%% INITIALIZATION %% 

debugFlag = 0; % for programmer debugging

% Resolve settings inputs; fall back to default values if necessary
if nargin<1, s = struct; end % create dummy structure if not present
if ~isfield(s,'re'), s.re = 'erbium'; end % default to erbium
if ~isfield(s,'arch'), s.arch = 'ring'; end % architecture can be (r)ing or (l)inear or (f)eedback ring
if ~isfield(s,'epsoutdb'), s.epsoutdb = 0; end % excess loss at output
switch lower(s.re(1))
    case 'e' % erbium
        if ~isfield(s,'L'), s.L = 10; end
        if ~isfield(s,'Ps'), s.Ps = 30e-6; end
        if ~isfield(s,'Pp'), s.Pp = 100e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1535e-9; end
        if ~isfield(s,'lamP'), s.lamP = 980e-9; end
        if ~isfield(s,'alpha'), s.alph = 6.5; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); end
        if ~isfield(s,'dCore'), s.dCore = 5.5e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end   
        tau32 = 10e-6;
        tau2 = 10e-3;
        [saP seP] = GetErSpectrum(s.lamP);
        [saS seS] = GetErSpectrum(s.lamS);
    case 'y' % ytterbium
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'Ps'), s.Ps = 10e-3; end
        if ~isfield(s,'Pp'), s.Pp = 250e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1064e-9; end
        if ~isfield(s,'lamP'), s.lamP = 975e-9; end
        if ~isfield(s,'alpha'), s.alph = 350; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(975)); end
        if ~isfield(s,'dCore'), s.dCore = 6e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end 
        tau32 = 1e-9;
        tau2 = 0.77e-3;   
        [saP seP] = GetYbSpectrum(s.lamP);
        [saS seS] = GetYbSpectrum(s.lamS);        
end
if ~isfield(s,'eps1db'), s.eps1db = 3; end
if ~isfield(s,'eps2db'), s.eps2db = 3; end
% fprintf('Running LaserPerformance with s.R2 = %d, s.R1 = %d, s.R = %d\n',round(s.R2*100),round(s.R1*100),round(s.R*100));
s = ResolveReflectances(s); % silly function to figure out which reflectances we got
            
%% DEFINE PARAMETERS %% 

% Rate equation parameters (Barnard Table I)
tau3 = tau32 + tau2;
h = 6.62606957e-34;
c = 3e8;
nuP = c/s.lamP;     % pump photon frequency
nuS = c/s.lamS;     % signal photon frequency

% Physical parameters (Barnard Table IIa)
b = s.dCore/2;      % radius of dopant in core
w = s.dField/2;     % radius of flux field travelling along core
Gamma = 1 - exp(-2*((b/w)^2));    % [1] overlap integral assuming Gaussian flux distribution
Ad = pi*b*b;                     % [m2] effective dopant area
etaP = tau3 / tau32;             % [1] pump efficiency
    % ** BELOW: MODIFIED FROM TABLE!! **
etaP = 1;
Beta2e = 1/(1 + tau32/tau2);     % [1]
Beta3e = Beta2e * (tau3/tau2);  % [1] 
    % modified from tau_32/tau_2 to tau_3/tau_2 so that Beta_3e --> 1 as
    % the pump level approaches the metastable level (see discussion 
    % in Barnard Section III-D)
etaQ = s.lamP / s.lamS;        % [1]

% Saturation paramters - Table III (three-level)
GammaP = Gamma;
GammaS = Gamma;
    % ** BELOW: CONCERNED ABOUT THESE VALUES
    % Strangely, these all seem to be modified by (1/h*nu) from the true
    % pump value ... is this intentional? all formulae using pump
    % saturation seem to consist of ratios so maybe it's okay ... can't
    % believe it should be like this though. I've modified to give real
    % pump values in Watts
    % Also concerned about cross-saturation values which 1) should be
    % symmetrical between PsCS and PpCS and 2) seems like the
    % cross-saturation should involve SOMETHING from the other "crossing" 
    % wavelength (not just all from the saturated wavelength)
PpIS = h*nuP*(Ad/GammaP) * (Beta2e/(etaP*tau2)) * (1/(saP+Beta3e*seP));
PpCS = h*nuP*(Ad/GammaP) * (1/tau2) * (1/saP);
PsIS = h*nuS*(Ad/GammaS) * (1/tau2) * (1/(saS+seS));
PsCS = h*nuS*(Ad/GammaS) * (Beta2e/(etaP*tau2)) * (1/(saS+Beta2e*seS));

% Laser parameters - Table IIb
alphaP = GammaP * s.Nt * saP;
alphaS = GammaS * s.Nt * saS;
delta = PsIS / PpIS;    % this is VERY sensitive to differences in 
                            % spectroscopy! in the end it doesn't matter 
                            % much as Gmax-->Inf for high N and long L
eps1 = 10^(-s.eps1db/10); % tx between output coupler and input of active fiber
eps2 = 10^(-s.eps2db/10); % tx between output of active fiber and the output coupler
switch lower(s.arch(1))
    case 'r' % ring laser
        kappa = s.R2;
    case 'l' % linear laser
        R1 = 1;
        R2 = s.R;
        R = sqrt(R1*R2);
    case 'f' % feedback ring
        kappa = s.R2;
        eps1 = eps1 * (1-kappa); 
end
eps = eps1 * eps2;

%% CALCULATE LASER PERFORMANCE %%

% Calculate total cavity losses (parasitic + output)
Gmax = exp(s.L*(alphaP/delta - alphaS)); % [ ]     
switch lower(s.arch(1))
    case 'l' % linear resonator      
        % can this be right? can result in value >1
        Tf = (1-eps2*eps2*R2) + (1-eps1*eps1*R1)*(eps2*eps2*R2)/(eps*R);
    case {'r','f'} % ring resonator
        Tf = (1-kappa*eps);
end

% Compute slope efficiency
switch lower(s.arch(1))
    case 'l' % linear resonator
        eta = etaQ * eps2 * ((1-R2)/(Tf)) * (PsIS / PsCS) * ...
            (1 - (Gmax * eps * R)^(-delta)); % from Table IV        
    case 'r' % ring resonator
        eta = etaQ * eps2 * ((1-kappa)/(Tf)) * (PsIS / PsCS) * ...
            (1 - (Gmax * eps * kappa)^(-delta)); % from Table IV
end
eta = eta * (10^(-s.epsoutdb/10)); % excess loss at output

% Compute lasing threshold
switch lower(s.arch(1))
    case 'l' % linear resonator
        Ppth = (1/etaQ) * PsCS * (alphaS * s.L - log(eps * R)) / ...
            (1 - (Gmax * eps * R)^(-delta)); % from Table IV        
    case 'r' % ring resonator
        Ppth = (1/etaQ) * PsCS * (alphaS * s.L - log(eps * kappa)) / ...
            (1 - (Gmax * eps * kappa)^(-delta)); % from Table IV 
end

% Compute round trip gain at threshold
%   Gain at lasing threshold must balance round trip losses, i.e. output
%   coupler losses and cavity losses
switch lower(s.arch(1))
    case 'l' % linear resonator
        Grt = 1 / (eps * R); % from discussion before eqn (25)
    case 'r' % ring resonator        
        Grt = 1 / (eps * kappa); % from discussion before eqn (25)
end

%% DISPLAY %%

if debugFlag
    % Print to prompt
    switch lower(s.arch(1))
        case 'l' % linear resonator
            fprintf('\nBarnard model laser type "%s": \n\t\t(R1 = %d%%, R2 = %d%%, eps_1 = %d%%, eps_2 = %d%%)\n\t\t eta = %.1f%% \t Pth = %.1f mW\n\n', ...
                scenario,round(100*R_1),round(100*R_2),round(100*eps_1),round(100*eps_2),100*eta,1000*P_p_th);
        case 'r' % ring resonator
            fprintf('\nBarnard model laser type "%s": \n\t\t(kappa = %d%%, eps_1 = %d%%, eps_2 = %d%%)\n\t\t eta = %.1f%% \t Pth = %.1f mW\n\n', ...
                scenario,round(100*kappa),round(100*eps_1),round(100*eps_2),100*eta,1000*P_p_th);
    end
end

% Check small-signal solution (neglecting ASE)
if debugFlag
    Ppin = Ppth;
    Ppout = 0;
    Psin = 1e-6;
    Psout = Grt*Psin;
    resP = Ppout - Ppin*exp(-alphaP*s.L + (Ppin-Ppout)/PpIS + (Psin-Psout)/PpCS);
    resS = Psout - Psin*exp(-alphaS*s.L + (Psin-Psout)/PsIS + (Ppin-Ppout)/PsCS);
end

% Correct outputs if necessary (indicates a non-lasing scenario)
if eta<0,
    eta = 0; 
    P_p_th = 0;
    G_rt = 0;
end

end

%%% END OF MAIN FUNCTION  %%%

function s = ResolveReflectances(s)
% figure out what the reflectances should be based on some user input

popupFlag = false;
if ~isfield(s,'R'), % R not given
    switch lower(s.arch(1))
        case 'l'
            popupFlag = true;
            if isfield(s,'R2'), % R2 given 
                if ~isfield(s,'R1'), s.R1 = 1; end % assume R1 is 1
                s.R = sqrt(s.R2*s.R1);
                popupFlag = false;
            elseif isfield(s,'R1') % R1 given but not R2
                s.R2 = s.R1; % assume the significant number is the output coupler
                s.R1 = 1;
                s.R = sqrt(s.R2*s.R1);
            else % neither R1 nor R2 given
                s.R2 = 0.5;
                s.R1 = 1;
                s.R = sqrt(s.R2*s.R1);
            end            
        case {'f','r'}
            if isfield(s,'R2'), s.R2 = s.R2;
            elseif isfield(s,'R1'), s.R2 = s.R1;
            else s.R2 = 0.5;
            end    
    end
elseif isfield(s,'R2') % R and R2 given
    switch lower(s.arch(1))
        case 'l'
            s.R1 = s.R*s.R/s.R2;  
            popupFlag = false;
        case {'f','r'}
            % ignore R; we'll just use R2
            s.R2 = s.R2;
    end
elseif isfield(s,'R1') % R and R1 given
    switch lower(s.arch(1))
        case 'l'
            s.R2 = s.R*s.R/s.R1;
            popupFlag = true;
        case {'f','r'}
            % ignore R; we'll just use R1
            s.R2 = s.R1;
    end     
    popupFlag = true;
else % only R given
    switch lower(s.arch(1))
        case 'l'
            s.R1 = 1;
            s.R2 = s.R;
            s.R = sqrt(s.R1*s.R2);
        case {'f','r'}
            s.R2 = s.R;
    end    
end

if popupFlag
    refl = inputdlg({'Please confirm: Output coupler reflectance (%):';'Back reflector reflectance (%):'}, ...
        'Reflectances',1,{num2str(100*s.R2),num2str(100*s.R1)});
    s.R2 = str2double(refl{1})/100; 
    s.R1 = str2double(refl{2})/100; 
    popupFlag = false;
end

end

