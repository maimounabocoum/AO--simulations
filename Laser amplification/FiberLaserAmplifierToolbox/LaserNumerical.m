function [Psout,hFig] = LaserNumerical(s)
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

if ~isfield(s,'plotFlag'), s.plotFlag = 1; end
plotFlag = s.plotFlag; 

if ~isfield(s,'re'), s.re = 'y'; end % default to erbium
if ~isfield(s,'arch'), s.arch = 'ring'; end % architecture can be (r)ing or (l)inear or (f)eedback ring
switch lower(s.re(1))
    case 'e' % erbium
        if ~isfield(s,'L'), s.L = 10; end
        if ~isfield(s,'Pp'), s.Pp = 100e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1535e-9; end
        if ~isfield(s,'lamP'), s.lamP = 980e-9; end
        if ~isfield(s,'alpha'), s.alph = 6.5; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); end
        if ~isfield(s,'dCore'), s.dCore = 5.5e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end   
    case 'y' % ytterbium
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'Pp'), s.Pp = 250e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1064e-9; end
        if ~isfield(s,'lamP'), s.lamP = 975e-9; end
        if ~isfield(s,'alpha'), s.alph = 350; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetYbSpectrum(975)); end
        if ~isfield(s,'dCore'), s.dCore = 6e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end 
end
if ~isfield(s,'eps1db'), s.eps1db = 3; end
if ~isfield(s,'eps2db'), s.eps2db = 3; end
% fprintf('Running LaserPerformance with s.R2 = %d, s.R1 = %d, s.R = %d\n',round(s.R2*100),round(s.R1*100),round(s.R*100));
s = ResolveReflectances(s); % silly function to figure out which reflectances we got
            
%% DEFINE PARAMETERS %% 

% Rate equation parameters (Barnard Table I)
h = 6.62606957e-34;
c = 3e8;
nuP = c/s.lamP;     % pump photon frequency
nuS = c/s.lamS;     % signal photon frequency

% Physical parameters (Barnard Table IIa)
b = s.dCore/2;      % radius of dopant in core
w = s.dField/2;     % radius of flux field travelling along core
Gamma = 1 - exp(-2*((b/w)^2));    % [1] overlap integral assuming Gaussian flux distribution
Ad = pi*b*b;                     % [m2] effective dopant area
etaQ = s.lamP / s.lamS;        % [1]

% Saturation paramters - Table III (three-level)
GammaP = Gamma;
GammaS = Gamma;

% Laser parameters - Table IIb
eps1 = 10^(-s.eps1db/10); % tx between output coupler and input of active fiber
eps2 = 10^(-s.eps2db/10); % tx between output of active fiber and the output coupler
switch lower(s.arch(1))
    case {'f','r'} % feedback ring, ring laser
        kappa = s.R2;
end
eps = eps1 * eps2;

%% CALCULATE LASER PERFORMANCE %%

% *** Define simulation parameters ***
maxIterationsAllowed = 25;           % [1] solution iterations allowed before abandoning attempt
maxChangeAllowed = 0.005;            % [dB] maximum allowable residual (error) at boundary conditions

% Simulate laser until convergence achieved
ii = 1;
N = 100;
maxDiff = 100;
PsCirc = 5e-3/eps1;
s.plotFlag = 0;
stressedFlag = 0;
while ii<=maxIterationsAllowed
    % simulate laser
    switch lower(s.arch(1))
        case {'r'} % ring
            % allow for losses between output coupler and gain medium
            PsCirc = PsCirc * eps1; % tx between output coupler and input of active fiber
            % simulate gain medium
            s.Ps = PsCirc;
            sUsed = s;
            if s.Ps<2e-3 && (s.Pp>100e-3 || L>0.5) % under "stressing" conditions, change the parameters a bit to get us up and running
                sUsed.Pp = min([s.Pp,100e-3+ii*10e-3]);
                sUsed.dz = 0.01; % very small resolution to ensure convergence
				stressedFlag = 1;
                [Psz,~,sUsed] = AmplifierPerformance(sUsed); % signal power as a function of z throughout fiber
            else
                [Psz,~,sUsed] = AmplifierPerformance(sUsed); % signal power as a function of z throughout fiber
				stressedFlag = 0;
            end
            if Psz == 0 % catch failure to converge
                jj = 1;
                while jj<=10
                    sUsed.dz = sUsed.dz/2;
                    [Psz,~,sUsed] = AmplifierPerformance(sUsed); % signal power as a function of z throughout fiber
                    if Psz>0, break; end
                    jj = jj + 1;
                end
            end
            PsCirc = Psz(end); % circulating power out of fiber
            % allow losses between gain medium and output coupler
            PsCirc = PsCirc * eps2; % tx between output of active fiber and the output coupler
            % output power from coupler
            PsOut(ii) = PsCirc * (1-kappa); % (1-kappa) leave the laser
            PsCirc = PsCirc * kappa;        % (kappa) remains in the laser and circulates
        case {'l'} % linear 
        case {'f'} % feedback ring
    end
    % print signal power
    fprintf('Output signal power is %.2f mW; %.2f mW reflected back into the cavity ... ',1000*PsOut(ii),1000*PsCirc);    
    % check convergence
    if ii>1 && ~stressedFlag
        maxDiffS = max(abs(10*log10(Psz./Psz0))); % maximum difference between this iteration's signal powers and the last one's
        fprintf('Max change is %.3f dB (<=%.3f dB for convergence).\n',maxDiffS,maxChangeAllowed);
        if maxDiffS<maxChangeAllowed, break; end        
    else
        fprintf('\n');
    end
    Psz0 = Psz;
    ii = ii + 1;
end % end simulation loop
s.plotFlag = plotFlag;

%% OUTPUT

% plot evolution of output power
if debugFlag
    figure;
    subplot(211); plot(1:length(PsOut),PsOut); title('Output Power vs Time');
    subplot(212); plot(1:length(Psz),Psz); title('Signal Power in Fiber');
end

% simulate powers in fiber
s.output = 'z'; % call for powers as a function of z
data = AmplifierPerformance(s); % signal power as a function of z throughout fiber
oldTitle = get(get(gca,'title'),'string');
newTitle = strrep(oldTitle,'DFA','DFL');
set(get(gca,'title'),'string',newTitle); % change plot to say "fiber laser" instead of "fiber amplifier"

% output values
z = data{1};
Psz = data{2};
Psout = PsOut(end);
hFig = gcf;

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

