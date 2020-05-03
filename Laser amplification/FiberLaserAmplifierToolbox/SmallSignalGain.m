function [data hFig] = SmallSignalGain(s)
% Analytical solution to small-signal single-pass gain in rare-earth 
% doped fiber amplifier
% 
% For saturated (i.e. large-signal) gain, run numerical solution in
% EDFASinglePassGain_Numerical
% 
% Following Payne et al, "Ytterbium-Doped Silica Fiber Lasers: Versatile
%   Sources for the 1-1.2 um Region", IEEE J. S. T. in Quantum Electronics, 
%   Vol 1. No 1.  1995
% 
% Also see:
%   (for validation)
%      Barnard et al, "Rare Earth Doped Fiber Amplifiers and
%       Lasers", IEEE J. Q. Electronics 1994,vol 30, no8, pp 1817-1830
%   (for parameter calculation, Erbium spectroscopy)
%      RP-Photonics.com articles: "Saturation Power", "Doping
%       Concentration", "Erbium-Doped Gain Media"
%   (for Ytterbium spectroscopy)
%      Marciante and Zuegel, "High-gain, polarization-preserving, 
%       Yb-doped fiber amplifier for low-duty-cycle pulse amplification." 
%       Applied optics 45.26 (2006): 6798-6804.
%
% Release date 17 May 2013, Luke Rumbaugh, Clarkson University
% 
% Inputs:
%   s -     settings structure as below:

%
% Outputs: 
%   G -     single pass gain in dB 
%           (nW x nP x nL) where nW is number of wavelengths, nP is number
%           of pump powers, and nL is number of fiber lengths to calculate

%% INITIALIZATION %% 

debugFlag = 0; % for programmer debugging

% Resolve settings inputs; fall back to default values if necessary
if nargin<1, s = struct; end % create dummy structure if not present
if ~isfield(s,'re'), s.re = 'erbium'; end % default to erbium
if ~isfield(s,'mode'), s.mode = 'wavelength'; end % plot vs wavelength or pump power
if ~isfield(s,'plotFlag'), s.plotFlag = 1; end
switch lower(s.re(1))
    case 'e' % erbium
        if ~isfield(s,'L'), 
            s.L = 10;
            multLflag = true;
        else
            multLflag = false;
        end 
        if ~isfield(s,'Ps'), s.Ps = 30e-6; end
        if ~isfield(s,'Pp'), 
            s.Pp = 100e-3; 
            multPflag = true;
        else
            multPflag = false;
        end
        if ~isfield(s,'lamS'), s.lamS = 1535e-9; end
        if ~isfield(s,'lamP'), s.lamP = 980e-9; end
        if ~isfield(s,'alpha'), s.alph = 6.5; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); end
        if ~isfield(s,'dCore'), s.dCore = 5.5e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end   
        wl = 1e-9 * (1450:5:1600);
    case 'y' % ytterbium
        if ~isfield(s,'L'), 
            s.L = 1;
            multLflag = true;
        else
            multLflag = false;
        end        
        if ~isfield(s,'Ps'), s.Ps = 10e-3; end
        if ~isfield(s,'Pp'), 
            s.Pp = 250e-3;
            multPflag = true;
        else
            multPflag = false;
        end
        if ~isfield(s,'lamS'), s.lamS = 1064e-9; end
        if ~isfield(s,'lamP'), s.lamP = 975e-9; end
        if ~isfield(s,'alpha'), s.alph = 350; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(975)); end
        if ~isfield(s,'dCore'), s.dCore = 6e-6; end
        if ~isfield(s,'Gamma'), s.Gamma = 0.722; end
        if ~isfield(s,'dField'), s.dField = s.dCore / s.Gamma; end 
        wl = 1e-9 * (1020:5:1100);
end

% Verify number of wavelengths (this affects calculation and display mode)
if isfield(s,'multPflag'), multPflag = s.multPflag; end
if isfield(s,'multLflag'), multLflag = s.multLflag; end
switch lower(s.mode(1))
    case 'p'
        nW = 1;
        wl = s.lamS;
        wlp = s.lamP;
        Pp = s.Pp * (1:1:100)/100;
        if multLflag
            L = s.L * (50:25:125)/100;
        else
            L = s.L;
        end
    otherwise
        nW = length(wl);
        wlp = s.lamP;
        if multPflag
            Pp = s.Pp * (50:10:120)/100;
        else
            Pp = s.Pp;
        end
        if multLflag
            L = s.L * (50:25:125)/100;
        else
            L = s.L;
        end    
end

% System parameters (user usually won't have to adjust)
A = pi*((s.dCore/2)^2); % m^2, doped core diameter
h = 6.62606957e-034;    % planck's number
c = 3e8;                % speed of light in vaccuum (not modified by index 
                        % of refraction for photon frequency)
nup = c/wlp;            % Hz, frequency of pump photons
switch lower(s.re(1))
    case 'y'
        [sal,sel] = GetYbSpectrum(wl);
        [sap,sep] = GetYbSpectrum(wlp);
        qe = 1; % 1, quantum efficiency of pump (1 for ytterbium, somewhat lower for erbium)
        tau = 0.77e-3; % s, spontaneous decay rate of excited ions
        fiberStr = 'Yb-Doped';
    case 'e'
        [sal,sel] = GetErSpectrum(wl);
        [sap,sep] = GetErSpectrum(wlp);
        tau = 8e-3; % s, spontaneous decay rate of excited ions
        qe = 0.85; % 1, quantum efficiency of pump (1 for ytterbium, somewhat lower for erbium)
        fiberStr = 'Er-Doped';
end

%% CALCULATION %% 

% Compute saturation power
PpSat = (h*nup*A/s.Gamma)/((sap+sep)*tau*qe); % saturation power for pump

% Calculate single-pass gain at each wavelength
%   For each fiber length
for ll = 1:numel(L)
    % For each pump power
    for pp = 1:numel(Pp)        
        % 1. Calculate power absorbed
        Ppl(pp) = SolveForPpz(Pp(pp),PpSat,-s.Gamma*s.Nt*sap*L(ll));
        Pa(pp) = Pp(pp)-Ppl(pp);
        % 2. Calculate gain for each wavelength (in dB)
        for ww = 1:numel(wl)
            G(ww,pp,ll) = (10*log10(exp(1)))*((qe*(sal(ww)+sel(ww))*tau*Pa(pp))./(h*nup*A/s.Gamma) - s.Gamma*s.Nt*sal(ww)*L(ll)); % dB, gain at each wavelength, at each pump power
        end
        % 3. Define legend entry
        legStr{pp} = sprintf('P_p = %d mW',round(1000*Pp(pp)));
    end
end

% Generate plots
if s.plotFlag
    switch lower(s.mode(1))
        case 'p'
            % Plot for a single-wavelength scenario
            %   Generates a single plot: Each length will be a different trace vs
            %   pump power
            cmap = hsv(length(L));
            hFig = figure('windowstyle','docked');
            for ll = 1:length(L);
                plot(1000*Pp,G(1,:,ll),'linewidth',2,'color',cmap(ll,:),'marker','o'); hold on;
                legStr{ll} = sprintf('%.2f m',L(ll));
            end
            hT = title(sprintf('Single-Pass Gain at %d nm\nUsing %s Fiber',round(wl*1e9),fiberStr)); set(hT,'fontsize',16);
            set(gca,'fontsize',14,'linewidth',2); grid on;
            xlabel(sprintf('Pump power (mW) at %d nm',round(wlp*1e9))); ylabel('Single-Pass Gain (dB)');
            legend(legStr,'Location','Best');        
        otherwise
            % Plot for a multiple-wavelength scenario
            %   Each fiber length is a different subplot
            %   Each pump power is a different trace vs wavelength on each subplot
            % Format figures and axes
            hFig = figure('windowstyle','docked');
            nAxX = max([1,floor(length(L)/2)]);      % number of rows of subplots
            nAxY = ceil(length(L)/nAxX);    % number of columns of subplots
            % Generate plots (one trace per pump power)
            for ll = 1:length(L)
                subplot(nAxX,nAxY,ll);
                plot(wl*1e9,G(:,:,ll),'linewidth',2);
                hT = title(sprintf('Active Fiber Length = %.2f m',L(ll))); set(hT,'fontsize',16);
                set(gca,'fontsize',14,'linewidth',2); grid on;
                xlabel ('Wavelength (nm)'); ylabel('Gain (dB)');
                hold on;
                ylim([0 1.2*max(max(G(:,:,ll)))]);
            end
            legend(legStr,'Location','EastOutside');
    end
    
    
end

% output
data = {G wl};

end
%%% END OF MAIN FUNCTION %%%


% SolveForPpz
% Solve for power at position L on fiber
function [solution,residual] = SolveForPpz(Pp0,Ps,k)
% based on equation: log(Pp(z) / Pp(0)) + (Pp(z) - Pp(0))/Ps = -N*sig_ap*a
% where k is the right hand side (should be negative)

if nargin<1
    Pp0 = 0.01;
    Ps = 0.002; % about right for 915 nm
    N = 1.5e25; % 550 ppm doping
    sap = 8e-25; % about right for 915 nm
    z = 1;
    k = -N*sap*z;
end

% initial guess is that between 0.01 and 0.99 of the power is absorbed
Ppz = Pp0*exp(k);
res = log(Ppz/Pp0) + (Ppz-Pp0)/Ps - k;
% solution = Ppz;
% residual = res;

% solving numerically, we can incorporate the saturation term as well
absDb = [0:0.001:60];
txLin = 10.^(-absDb./10);
Ppz = txLin*Pp0; 

% calculate residual for each Ppz guess
res = log(Ppz./Pp0) + (Ppz-Pp0)./Ps - k;

% output value with lowest residual
[minres,iimin] = min(abs(res));
solution = Ppz(iimin);
residual = minres;

% debugging
%   plot(1000*(Pp0-Ppz),res); xlabel('Power absorbed in fiber (mW)'); ylabel('Residual error in calc'); title(sprintf('Absorbed power in fiber launching %d mW pump',1000*Pp0));
%   fprintf('Pump power %f mW \t Term 1 %f \t Term 2 %f \t Term 3 %f \t Residual %f \t Pump at z %f mW \t Absorbed power %f mW\n',Pp0,log(Ppz(iimin)./Pp0),(Ppz(iimin)-Pp0)./Ps,-k,residual,solution,Pp0-solution);

end
