function [data,hFig,s] = AmplifierPerformance(s)
% [data,hFig,s] = AmplifierPerformance(s)
% Numerical solution to the powers propagating in a YDFA
% Using finite difference methods to solve the coupled PDEs
% 
% s = struct('re','y','pump','core','Ps',0,'Pp',50e-3,'dz',0.01,'L',0.5);
%
% Uses Forward-Time, Centered-Space explicit finite difference method
% 
% References: 
%   "Strongly pumped fiber lasers", Kelson & Hardy, 1998
%   "Relaxation methods for partial differential equations: Applications to
%   electrostatics", David Robertson, 2010
%   OASiX v. 4.0 (validation tool)
%
% Revisions 
%   Based on 25 April 13 version of ODE solver
%   Created 31 May 13
% 
% Inputs:
%   lamP - wavelength of pump
%   lamS - wavelength of signal
%   Pp - constant power of pump input
%   Ps - constant power of signal input
%   L - amplifier length
%   pump - pumping scheme: "core" or "cladding" options expected
%   s - settings structure with fields:
%       direction - pump direction: forward ("f"), reverse ("r"),
%           bidrectional ("b"), or double-pass ("2") options expected
%       R - length adjustment factor (for faster simulation)
%       diaP - diameter of fiber where pump travels (core or inner
%           cladding, depending on pumping scheme)
%       diaS - diameter of fiber where signal travels (doped core)
%       Ntot - doping concentration
%       dlam - optical sampling bandwidth (for ASE)
%       displayFlag - generate plots?
%       dz - spatial resolution of calculation
%
% Outputs:
%   Z - vector of positions along fiber, for which solution was calculated
%   Y - matrix of values along fiber, with columns [P,S,ASEf,ASEb,n2]:
%       P - power of pump (direction determined from s.direction)
%       S - power of signal (forward propagating power at lamS)
%       ASEf - spectrally integrated forward propagating ASE (will not include lamS)
%       ASEb - spectrally integrated backward propagating ASE (will include lamS)
%       n2 - inversion RATIO (i.e. domain = [0,1]) where n2 = N2/Ntot

debugFlag = 0;

%% INITIALIZATION %%

% Resolve settings inputs; fall back to default values if necessary
if nargin<1, s = struct; end % create dummy structure if not present
if ~isfield(s,'re'), s.re = 'erbium'; end % default to erbium
if ~isfield(s,'pump'), s.pump = 'core'; end % default to core pumped
if ~isfield(s,'plotFlag'),s.plotFlag = 1; end     % [bool] generate plots at end?
if ~isfield(s,'mode'), s.mode = 'power'; end % can be power, gain, or ASE
if ~isfield(s,'direction'), s.direction = 'f';  end     % pumping direction
if ~isfield(s,'dlam'),s.dlam = 2e-9;            end     % [m] sampling bandwidth
if ~isfield(s,'Pase'), s.Pase = 0; end      % ASE seed power
switch [lower(s.re(1)) lower(s.pump)]
    case 'ecore' % erbium core-pumped
        if ~isfield(s,'L'), s.L = 10; end
        if ~isfield(s,'Ps'), s.Ps = 30e-6; end
        if ~isfield(s,'Pp'), s.Pp = 100e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1535e-9; end
        if ~isfield(s,'lamP'), s.lamP = 980e-9; end
        if ~isfield(s,'alpha'), s.alph = 6.5; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); end
        if ~isfield(s,'dCore'), s.dCore = 5.5e-6; end
        if ~isfield(s,'tau'),s.tau = 10e-3;              end     % [ms] upper state lifetime    
		if ~isfield(s,'nP'),
			s.nP = min([75,s.L/0.01]); % number of points on line: default to 100, but don't go below 1 cm resolution
			s.nP = max([s.nP,s.L/0.10]); % ... and, don't go above 10 cm resolution
		end
        lam = [s.lamP (1450e-9:s.dlam:1600e-9)]';% [m] wavelength vector
        [sig12 sig21] = GetErSpectrum(lam);          
        ampName = 'EDFA';
    case 'eclad' % erbium cladding-pumped
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'Ps'), s.Ps = 30e-6; end
        if ~isfield(s,'Pp'), s.Pp = 100e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1535e-9; end
        if ~isfield(s,'lamP'), s.lamP = 980e-9; end
        if ~isfield(s,'alpha'), s.alph = 6.5; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetErSpectrum(1530)); end
        if ~isfield(s,'dCore'), s.dCore = 5.5e-6; end
        if ~isfield(s,'tau'),s.tau = 10e-3;              end     % [ms] upper state lifetime    
		if ~isfield(s,'nP'),
			s.nP = min([75,s.L/0.01]); % number of points on line: default to 100, but don't go below 1 cm resolution
			if s.Pp>10
				s.nP = max([s.nP,s.L/0.02]); % ... and, don't go above 2 cm resolution for high power pumping
			else
				s.nP = max([s.nP,s.L/0.1]); % or 10 cm resolution for low power pumping
			end
		end		
        lam = [s.lamP (1450e-9:s.dlam:1600e-9)]';% [m] wavelength vector
        [sig12 sig21] = GetErSpectrum(lam);          
        ampName = 'EDFA';        
    case 'ycore' % ytterbium
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'Ps'), s.Ps = 10e-3; end
        if ~isfield(s,'Pp'), s.Pp = 250e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1064e-9; end
        if ~isfield(s,'lamP'), s.lamP = 975e-9; end
        if ~isfield(s,'alpha'), s.alph = 350; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetYbSpectrum(975)); end
        if ~isfield(s,'dCore'), s.dCore = 6e-6; end
        if ~isfield(s,'tau'),s.tau = 0.77e-3;              end     % [ms] upper state lifetime
		if ~isfield(s,'nP'),
			s.nP = min([75,s.L/0.01]); % number of points on line: default to 100, but don't go below 1 cm resolution
			s.nP = max([s.nP,s.L/0.10]); % ... and, don't go above 10 cm resolution
		end		
        lam = [s.lamP (1020e-9:s.dlam:1100e-9)]';% [m] wavelength vector
        [sig12 sig21] = GetYbSpectrum(lam);   
        ampName = 'YDFA';
    case 'yclad' % ytterbium
        if ~isfield(s,'L'), s.L = 1; end
        if ~isfield(s,'Ps'), s.Ps = 10e-3; end
        if ~isfield(s,'Pp'), s.Pp = 250e-3; end
        if ~isfield(s,'lamS'), s.lamS = 1064e-9; end
        if ~isfield(s,'lamP'), s.lamP = 975e-9; end
        if ~isfield(s,'alpha'), s.alph = 350; end
        if ~isfield(s,'Nt'), s.Nt = ConvAbsDB2N(s.alph,GetYbSpectrum(975)); end
        if ~isfield(s,'dCore'), s.dCore = 6e-6; end
        if ~isfield(s,'tau'),s.tau = 0.77e-3;              end     % [ms] upper state lifetime
		if ~isfield(s,'nP'),
			s.nP = min([75,s.L/0.01]); % number of points on line: default to 100, but don't go below 1 cm resolution
			if s.Pp>10
				s.nP = max([s.nP,s.L/0.02]); % ... and, don't go above 2 cm resolution for high power pumping
			else
				s.nP = max([s.nP,s.L/0.1]); % or 10 cm resolution for low power pumping
			end
		end			
        lam = [s.lamP (1020e-9:s.dlam:1100e-9)]';% [m] wavelength vector
        [sig12 sig21] = GetYbSpectrum(lam);   
        ampName = 'YDFA';     
end
if ~isfield(s,'dz'),s.dz = s.L/s.nP;                 end     % [m] positional solution mesh resolution
s.nP = s.L/s.dz; % recalculate number of points
%   resolve overlap parameters
if ~isfield(s,'Gamma'), s.Gamma = 0.722; end % assumes 80% core overlap with Gaussian distributed intensity cross-section
switch lower(s.pump)
    case 'core'
        s.GammaP = s.Gamma;                 % [1] effective spatial mode overlap of pump over doped core
        s.GammaS = s.Gamma;                             % [1] effective spatial mode overlap of signal over doped core
        s.dFieldP = s.dCore / s.GammaP;                 
        s.dFieldS = s.dCore / s.GammaS;
        Ap = pi*(s.dCore/2)^2;               % [m2] cladding area of pump waveguide
        As = pi*(s.dCore/2)^2;               % [m2] core area of signal waveguide
    case 'clad'
        if ~isfield(s,'dClad'), s.dClad = 125e-6; end
        Ap = pi*(s.dClad/2)^2;               % [m2] cladding area of pump waveguide
        As = pi*(s.dCore/2)^2;               % [m2] core area of signal waveguide
        s.GammaP = s.Gamma * (As/Ap);         % [1] effective spatial mode overlap of pump over doped core
        s.GammaS = s.Gamma;                             % [1] effective spatial mode overlap of signal over doped core
        s.dFieldP = s.dClad / s.GammaP;                 
        s.dFieldS = s.dCore / s.GammaS;
end

% *** Define simulation parameters ***
if ~isfield(s,'diffMethod'), s.diffMethod = 'c'; end % central difference
maxIterationsAllowed = 1e3;           % [1] solution iterations allowed before abandoning attempt
makeInitialGuessFlag = false;        % [bool] perform simple calculation to initialize fiber with a reasonable solution guess; this gives >5x speedup
maxChangeAllowed = 0.005;            % [dB] maximum allowable residual (error) at boundary conditions
if debugFlag, figure; end

    % define physical constants
h = 6.626e-34;                      % [J-s] Planck's constant
c = 3e8;                            % [m/s] speed of light vacuum
n = 1.5;                            % [1] index of refraction for ZBLAN fiber
vg = c/n;
    % define system parameters
z = 0:s.dz:s.L;                         % [m] positional solution mesh
dt = s.dz/(2*vg);
lambdaF = dt/s.dz;                  % for forward propagating signals
lambdaB = dt/(-s.dz);               % for backward propagating signals
dzSignF = 1;                        % for forward propagating signals
dzSignB = -1;                       % for backward propagating signals
alpha = 3e-3;                       % [1/m] attenuation coefficient
    % derived parameters
v = c./lam;                          % [1/s] frequency corresponding to wavelengths
iiP = 1;                            % [#] index of pump wavelength in lam vector
[~,iiS] = min(abs(s.lamS-lam));       % [#] index of signal wavelength in lam vector
aseTerm = GetASETerm(s.GammaS,sig21,lam,h,c,s.dlam,s.Nt); % term to handle ASE; function of wavelength alone
% define solution matrix
Nlam = length(lam);                 % number of wavelengths
Nz = length(z);                     % number of spatial mesh nodes
Y = zeros(Nz,2*Nlam);               % each row is a spatial node; each column is a wavelength. bidirectional propagation for all wavelengths including pump.
cl.Pf = 1;                          % column number of pump forward
cl.Pb = 2;                          % column number of pump backwards
cl.Sf = iiS+1;                      % column number of signal forward
cl.Sb = Nlam+iiS+1;                 % column number of signal backward (double-pass case only)
cl.ASEf = [3:iiS,(iiS+2):(Nlam+1)]; % column numbers of ASE forward
cl.ASEb = [(Nlam+2):(Nlam+iiS),(Nlam+iiS+2):(2*Nlam)]; % column numbers of ASE backward
iiASE = cl.ASEf - 1;                % index number for coefficients etc
% define convergence vectors
dPdt = zeros(1,Nz);
dSfdt = zeros(1,Nz);
dAfdt_max = zeros(1,Nz);
dAbdt_max = zeros(1,Nz);


%% Compute Solution

fprintf('Solving BVP ... '); tic;   % notify user and initialize clock
fprintf('(L=%dcm,Pp=%dmW,Ps=%duW)(#points=%d,dz=%.1fcm) ... ',round(100*s.L),round(1000*s.Pp),round(1e6*s.Ps),s.nP,s.dz*100); % nice to see details about simulation
maxres = 1;                         % track the maximum residual
n = 1;                              % count the iterations performed

% *** Apply initial conditions ***
N2 = s.Nt * zeros(Nz,1);
N1 = s.Nt * ones(Nz,1);

% *** Apply boundary conditions *** 
%   Pump powers
switch s.direction
	case 'f'
		dzSignP = dzSignF;
        cl.P = cl.Pf;
        lambdaP = lambdaF;
        Y(1,cl.Pf) = s.Pp;        % pump power forward at z = 0
        Y(end,cl.Pb) = 0;       % pump power backward at z = L
	case {'b','r'}
		dzSignP = dzSignB;
        cl.P = cl.Pb;
        lambdaP = lambdaB;
        Y(1,cl.Pf) = 0;        % pump power forward at z = 0
        Y(end,cl.Pb) = s.Pp;       % pump power backward at z = L
end
%   Signal powers
Y(1,cl.Sf) = s.Ps;        % signal power forward at z = 0
Y(end,cl.Sb) = 0;       % signal power backward at z = L
%   ASE powers
Y(1,cl.ASEf) = s.Pase;       % ASE power forward at z = 0
Y(end,cl.ASEb) = s.Pase;     % ASE power backward at z = L

% *** Solve fiber *** 
%   Gauss-Siedel method
while n <= maxIterationsAllowed
    
    % For each update step, propagate boundary condition through the whole fiber
    %   0. Points are numbered as:
    %       1 : boundary condition at z = 0
    %       2 through Nz-1 : interior points
    %       Nz : boundary value at z = L
    %   1. Use Forward-Time Centered-Space stencil for points 2 thru Nz-1
    %       P_n+1,z = P_n,z + v*dt*(Alpha*P - (P_n,z+1 - P_n,z-1)/(2*dz))  
    %   2. Use Forward-Time Backward-Space stencil for point Nz
    %       P_n+1,z = P_n,z + v*dt*(Alpha*P - (P_n,z - P_n,z-1)/(dz))
        
    % Step through interior points
    for zz = 2:(Nz-1) % index zz is the number of the point to be updated
       
        zzF = zz;
        zzB = (Nz-(zzF-1));
        
        %%% Pump: FTCS %%%
        %   points to be used for explicit calculation
        switch s.direction(1)
            case 'f' % propagating forward from beginning to end
                ii0 = zzF - 1;                  % last updated point 
                ii1 = zzF;                      % this point 
                ii2 = zzF + 1;                  % next point to update
            case {'r','b'} % propagating backward from end to beginning             
                ii0 = zzB + 1;                  % last updated point 
                ii1 = zzB;                      % this point
                ii2 = zzB - 1;                  % next point to update   
        end
        P0 = Y(ii0,cl.P);    % previous point, recently updated w valid information
        P1 = Y(ii1,cl.P);    % this point, prior to update
        P2 = Y(ii2,cl.P);    % next point, prior to update
        switch lower(s.diffMethod(1))
            case 'c'
                %   inhomogenous term defining physics of fiber
                D0 = - s.GammaP*sig12(iiP)*N1(ii0)*P0 ...   % absorption term
                    + s.GammaP*sig21(iiP)*N2(ii0)*P0 ...   % emission term
                    - alpha*P0;                           % attenuation term
                D1 = - s.GammaP*sig12(iiP)*N1(ii1)*P1 ...   % absorption term					
                    + s.GammaP*sig21(iiP)*N2(ii1)*P1 ...   % emission term
                    - alpha*P1;                           % attenuation term
                D2 = - s.GammaP*sig12(iiP)*N1(ii2)*P2 ...   % absorption term
                    + s.GammaP*sig21(iiP)*N2(ii2)*P2 ...   % emission term
                    - alpha*P2;                           % attenuation term
                %   centered-space differential terms
                dPdz1 = dzSignP*(-vg*lambdaP/2)*(P2 - P0);
                dPdz2 = (vg*vg*lambdaP*lambdaP/2)*(P2 - 2*P1 + P0);
                %   differentials of the inhomogenous
                Dscaled = (vg*dt)*D1;
                dDdt = 0*(vg*dt/2)*(D1 + D1);
                dDdz = dzSignP*(-vg*vg*dt*lambdaP/4)*(D2 - D0);
                %   update by Taylor series expansion
                P1new = P1 + dPdz1 + dPdz2 + Dscaled + dDdt + dDdz;
                %   recalculate using updated diffusion term
%                 D1new = - s.GammaP*sig12(iiP)*N1(zz)*Pf1new ...   % absorption term
%                     + s.GammaP*sig21(iiP)*N2(zz)*Pf1new ...   % emission term
%                     - alpha*Pf1new;                           % attenuation term
%                 dDdt = (vg*dt/2)*(D1new - D1);
%                 Pf1new = Pf1 + dPfdz1 + dPfdz2 + Dscaled + dDdt + dDdz;
            case {'f','b'}                
                de = - s.GammaP*sig12(iiP)*N1(ii0)*P0 ...   % absorption term
                    + s.GammaP*sig21(iiP)*N2(ii0)*P0 ...   % emission term
                    - alpha*P0;                           % attenuation term
                dPfdz = (P1 - P0) / (s.dz);                
                %   update by forward-marching time differential
                P1new = P1 + vg*dt*(de - dPfdz);
        end        
        %   update the matrix
%         Pf1new = max([0,Pf1new]); % disallow negative values
        Y(ii1,cl.P) = P1new;
        %   note change at this point (checking convergence)
        dPdt(ii1) = abs(10*log10(P1new/P1));
        if isnan(dPdt(ii1)) || isinf(dPdt(ii1)), dPdt(ii1) = 0; end
        
        %%% Signal: FTCS %%%
        %   points to be used for explicit calculation
        Sf0 = Y(zzF-1,cl.Sf);    % previous point zz-1, recently updated
        Sf1 = Y(zzF,cl.Sf);      % this point, prior to update
        Sf2 = Y(zzF+1,cl.Sf);    % next point zz+1, prior to update
        switch lower(s.diffMethod(1))
            case 'c'
                %   inhomogenous term defining physics of fiber
                D0 = - s.GammaS*sig12(iiS)*N1(zzF-1)*Sf0 ...   % absorption term
                    + s.GammaS*sig21(iiS)*N2(zzF-1)*Sf0 ...   % stimulated emission term
                    + aseTerm(iiS)*N2(zzF-1)/s.Nt ...         % spontaneous emission term
                    - alpha*Sf0;                           % attenuation term
                D1 = - s.GammaS*sig12(iiS)*N1(zzF)*Sf1 ...   % absorption term
                    + s.GammaS*sig21(iiS)*N2(zzF)*Sf1 ...   % emission term
                    + aseTerm(iiS)*N2(zzF)/s.Nt ...         % spontaneous emission term
                    - alpha*Sf1;                           % attenuation term
                D2 = - s.GammaS*sig12(iiS)*N1(zzF+1)*Sf2 ...   % absorption term
                    + s.GammaS*sig21(iiS)*N2(zzF+1)*Sf2 ...   % emission term
                    + aseTerm(iiS)*N2(zzF+1)/s.Nt ...         % spontaneous emission term
                    - alpha*Sf2;                           % attenuation term
                %   centered-space differential terms
                dSfdz1 = dzSignF*(-vg*lambdaF/2)*(Sf2 - Sf0);
                dSfdz2 = (vg*vg*lambdaF*lambdaF/2)*(Sf2 - 2*Sf1 + Sf0);
                %   differentials of the inhomogenous
                Dscaled = (vg*dt)*D1;
                dDdt = 0*(vg*dt/2)*(D1 + D1);
                dDdz = dzSignF*(-vg*vg*dt*lambdaF/4)*(D2 - D0);
                %   update by Taylor series expansion
                Sf1new = Sf1 + dSfdz1 + dSfdz2 + Dscaled + dDdt + dDdz;
                %   recalculate using updated diffusion term
%                 D1new = - s.GammaS*sig12(iiS)*N1(zz)*Sf1new ...   % absorption term
%                     + s.GammaS*sig21(iiS)*N2(zz)*Sf1new ...   % emission term
%                     - alpha*Sf1new;                           % attenuation term
%                 dDdt = (vg*dt/2)*(D1new - D1);
%                 Sf1new = Sf1 + dSfdz1 + dSfdz2 + Dscaled + dDdt + dDdz;
%                 Sf1new = Sf1 + dSfdz1 + Dscaled;
            case {'f','b'}                
                de = - s.GammaS*sig12(iiS)*N1(zzF)*Sf0 ...   % absorption term
                    + s.GammaS*sig21(iiS)*N2(zzF)*Sf0 ...   % stimulated emission term
                    + aseTerm(iiS)*N2(zzF)/s.Nt ...         % spontaneous emission term
                    - alpha*Sf0;                           % attenuation term
                %   forward-space discretized spatial differential
                dSfdz = (Sf1 - Sf0) / (s.dz);
                %   update by forward-marching time differential
                Sf1new = Sf1 + vg*dt*(de - dSfdz);
        end        
        %   update the matrix
%         Sf1new = max([0,Sf1new]); % disallow negative values
        Y(zzF,cl.Sf) = Sf1new;
        %   note change at this point (checking convergence)
        dSfdt(zzF) = abs(10*log10(Sf1new/Sf1));
        if isnan(dSfdt(zzF)) || isinf(dSfdt(zzF)), dSfdt(zzF) = 0; end        
                        
        %%% ASE Forward: FTCS %%%
        %   points to be used for explicit calculation
        Af0 = Y(zzF-1,cl.ASEf);    % previous point zz-1, recently updated
        Af1 = Y(zzF,cl.ASEf);      % this point, prior to update
        Af2 = Y(zzF+1,cl.ASEf);    % next point zz+1, prior to update
        switch lower(s.diffMethod(1))
            case 'c'
                %   inhomogenous term defining physics of fiber
                D0 = - s.GammaS*sig12(iiASE)'.*N1(zzF-1).*Af0 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzF-1).*Af0 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzF-1)./s.Nt ...         % spontaneous emission term
                    - alpha*Af0;                           % attenuation term
                D1 = - s.GammaS*sig12(iiASE)'.*N1(zzF).*Af1 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzF).*Af1 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzF)./s.Nt ...         % spontaneous emission term
                    - alpha*Af1;                           % attenuation term
                D2 = - s.GammaS*sig12(iiASE)'.*N1(zzF+1).*Af2 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzF+1).*Af2 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzF+1)./s.Nt ...         % spontaneous emission term
                    - alpha*Af2;                           % attenuation term
                %   centered-space differential terms
                dAfdz1 = dzSignF*(-vg*lambdaF/2)*(Af2 - Af0);
                dAfdz2 = (vg*vg*lambdaF*lambdaF/2)*(Af2 - 2*Af1 + Af0);
                %   differentials of the inhomogenous terms
                Dscaled = (vg*dt)*D1;
                dDdt = 0*(vg*dt/2)*(D1 + D1);
                dDdz = dzSignF*(-vg*vg*dt*lambdaF/4)*(D2 - D0);
                %   update by Taylor series expansion
                Af1new = Af1 + dAfdz1 + dAfdz2 + Dscaled + dDdt + dDdz;
                %   recalculate using updated diffusion term
%                 D1new = - s.GammaS*sig12(iiS)*N1(zz)*Sf1new ...   % absorption term
%                     + s.GammaS*sig21(iiS)*N2(zz)*Sf1new ...   % emission term
%                     - alpha*Sf1new;                           % attenuation term
%                 dDdt = (vg*dt/2)*(D1new - D1);
%                 Sf1new = Sf1 + dSfdz1 + dSfdz2 + Dscaled + dDdt + dDdz;
%                 Sf1new = Sf1 + dSfdz1 + Dscaled;
            case {'f','b'}                
                de = -s.GammaS*sig12(iiASE)'.*N1(zzF).*Af0 ...
                    + s.GammaS*sig21(iiASE)'.*N2(zzF).*Af0 ...
                    + aseTerm(iiASE)'.*N2(zzF)./s.Nt ...
                    - alpha*Af0; % difference between steps
                %   forward-space discretized spatial differential
                dAfdz = (Af1 - Af0) ./ (s.dz);
                %   update by forward-marching time differential
                Af1new = Af1 + vg*dt*(de - dAfdz);
        end             
        Y(zzF,cl.ASEf) = Af1new; 
        dAfdt = abs(10*log10(Af1new./Af1));  
        dAfdt_max(zzF) = max(dAfdt);
        if isnan(dAfdt_max(zzF)) || isinf(dAfdt_max(zzF)), dAfdt_max(zzF) = 100; end       
        
        %%% ASE Backward: FTCS %%%
        %   points to be used for explicit calculation
        %   propagating from end to beginning        
        Ab0 = Y(zzB+1,cl.ASEb);    % previous point zz+1, recently updated
        Ab1 = Y(zzB,cl.ASEb);      % this point, prior to update
        Ab2 = Y(zzB-1,cl.ASEb);    % next point zz-1, prior to update
        switch lower(s.diffMethod(1))
            case 'c'
                %   inhomogenous term defining physics of fiber
                D0 = - s.GammaS*sig12(iiASE)'.*N1(zzB-1).*Ab0 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzB-1).*Ab0 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzB-1)./s.Nt ...         % spontaneous emission term
                    - alpha*Ab0;                           % attenuation term
                D1 = - s.GammaS*sig12(iiASE)'.*N1(zzB).*Ab1 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzB).*Ab1 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzB)./s.Nt ...         % spontaneous emission term
                    - alpha*Ab1;                           % attenuation term
                D2 = - s.GammaS*sig12(iiASE)'.*N1(zzB+1).*Ab2 ...   % absorption term
                    + s.GammaS*sig21(iiASE)'.*N2(zzB+1).*Ab2 ...   % stimulated emission term
                    + aseTerm(iiASE)'.*N2(zzB+1)./s.Nt ...         % spontaneous emission term
                    - alpha*Ab2;                           % attenuation term
                %   centered-space differential terms
                dAbdz1 = dzSignB*(-vg*lambdaB/2)*(Ab2 - Ab0); 
                dAbdz2 = (vg*vg*lambdaB*lambdaB/2)*(Ab2 - 2*Ab1 + Ab0);
                %   differentials of the inhomogenous terms
                Dscaled = (vg*dt)*D1;
                dDdt = 0*(vg*dt/2)*(D1 + D1);
                dDdz = dzSignB*(-vg*vg*dt*lambdaB/4)*(D2 - D0); 
                %   update by Taylor series expansion
                Ab1new = Ab1 + dAbdz1 + dAbdz2 + Dscaled + dDdt + dDdz;
                %   recalculate using updated diffusion term
%                 D1new = - s.GammaS*sig12(iiS)*N1(zz)*Sf1new ...   % absorption term
%                     + s.GammaS*sig21(iiS)*N2(zz)*Sf1new ...   % emission term
%                     - alpha*Sf1new;                           % attenuation term
%                 dDdt = (vg*dt/2)*(D1new - D1);
%                 Sf1new = Sf1 + dSfdz1 + dSfdz2 + Dscaled + dDdt + dDdz;
%                 Sf1new = Sf1 + dSfdz1 + Dscaled;
            case {'f','b'}                
                de = -s.GammaS*sig12(iiASE)'.*N1(zzB).*Ab0 ...
                    + s.GammaS*sig21(iiASE)'.*N2(zzB).*Ab0 ...
                    + aseTerm(iiASE)'.*N2(zzB)./s.Nt ...
                    - alpha*Ab0; % difference between steps
                %   forward-space discretized spatial differential
                dAbdz = (Ab1 - Ab0) ./ (s.dz);
                %   update by forward-marching time differential
                Ab1new = Ab1 + vg*dt*(de - dAbdz);
        end             
        Y(zzB,cl.ASEb) = Ab1new; 
        dAbdt = abs(10*log10(Ab1new./Ab1));  
        dAbdt_max(zzB) = max(dAbdt);
        if isnan(dAbdt_max(zzB)) || isinf(dAbdt_max(zzB)), dAbdt_max(zzB) = 100; end           
        
        % deal ASE out of convergence condition for now      
        dAfdt_max(zzF) = 0;
        dAbdt_max(zzB) = 0;
        
    end
    
    % Now compute the very last point using a [forward/backward] difference
    % (i.e., just using two spatial points, not three, since we're up 
    % against the boundary)
    for zz = Nz
        
        zzF = zz;
        zzB = (Nz-(zzF-1));
        
        %   points to be used for explicit calculation
        switch s.direction(1)
            case 'f' % propagating forward from beginning to end
                ii0 = zzF - 1;                  % last updated point 
                ii1 = zzF;                      % this point 
            case {'r','b'} % propagating backward from end to beginning             
                ii0 = zzB + 1;                  % last updated point 
                ii1 = zzB;                      % this point
        end                
        % Pump: FT[F/B]S
        %   points to be used for explicit calculation
        P0 = Y(ii0,cl.P);       % last updated point (valid data)
        P1 = Y(ii1,cl.P);       % this point prior to updating
        %   differential equation per physics of fiber
        de = - s.GammaP*sig12(iiP)*N1(zz)*P0 ...   % absorption term
            + s.GammaP*sig21(iiP)*N2(zz)*P0 ...   % emission term
            - alpha*P0;                           % attenuation term
        %   backward-space discretized spatial differential
        dPfdz = dzSignP*(P1 - P0)/(dzSignP*s.dz);
        %   update by forward-marching time differential
        P1new = P1 + vg*dt*(de - dPfdz);
        P1new = max([0,P1new]); % disallow negative values
        Y(ii1,cl.P) = P1new;
        %   note change at this point (checking convergence)
        dPdt(ii1) = abs(10*log10(P1new/P1));
        if isnan(dPdt(ii1)) || isinf(dPdt(ii1)), dPdt(ii1) = 0; end
        
        % Signal: FTFS
        %   points to be used for explicit calculation
        Sf0 = Y(zzF-1,cl.Sf);    % previous point zz-1, recently updated
        Sf1 = Y(zzF,cl.Sf);      % this point, prior to update
        %   differential equation per physics of fiber
        de = - s.GammaS*sig12(iiS)*N1(zzF)*Sf0 ...   % absorption term
            + s.GammaS*sig21(iiS)*N2(zzF)*Sf0 ...   % stimulated emission term
            + aseTerm(iiS)*N2(zzF)/s.Nt ...         % spontaneous emission term
            - alpha*Sf0;                           % attenuation term
        %   centered-space discretized spatial differential
        dSfdz = dzSignF*(Sf1 - Sf0)/(dzSignF*s.dz);
        %   update by forward-marching time differential
        Sf1new = Sf1 + vg*dt*(de - dSfdz);
        Sf1new = max([0,Sf1new]); % disallow negative values
        Y(zzF,cl.Sf) = Sf1new;
        %   note change at this point (checking convergence)
        dSfdt(zzF) = abs(10*log10(Sf1new/Sf1));
        if isnan(dSfdt(zzF)) || isinf(dSfdt(zzF)), dSfdt(zzF) = 0; end
        
        % ASE Forward: FTFS
        %   points to be used for explicit calculation
        Af0 = Y(zzF-1,cl.ASEf);    % previous point zz-1, recently updated
        Af1 = Y(zzF,cl.ASEf);      % this point, prior to update
        %   differential equation per physics of fiber
        de = - s.GammaS.*sig12(iiASE)'.*N1(zzF).*Af0 ...   % absorption term
            + s.GammaS.*sig21(iiASE)'.*N2(zzF).*Af0 ...   % stimulated emission term
            + aseTerm(iiASE)'.*N2(zzF)./s.Nt ...         % spontaneous emission term
            - alpha*Af0;                           % attenuation term
        %   centered-space discretized spatial differential
        dAfdz = dzSignF*(Af1 - Af0)/(dzSignF*s.dz);
        %   update by forward-marching time differential
        Af1new = Af1 + vg*dt*(de - dAfdz);
        Af1new(Af1new<0) = 0; % disallow negative values
        Y(zzF,cl.ASEf) = Af1new;
        %   note change at this point (checking convergence)
        dAfdt = abs(10*log10(Af1new./Af1));  
        dAfdt_max(zzF) = max(dAfdt);
        if isnan(dAfdt_max(zzF)) || isinf(dAfdt_max(zzF)), dAfdt_max(zzF) = 100; end       
         
        % ASE Backward: FTBS
        %   points to be used for explicit calculation
        Ab0 = Y(zzB+1,cl.ASEb);    % previous point zz-1, recently updated
        Ab1 = Y(zzB,cl.ASEb);      % this point, prior to update
        %   differential equation per physics of fiber
        de = - s.GammaS*sig12(iiASE)'.*N1(zzB).*Ab0 ...   % absorption term
            + s.GammaS*sig21(iiASE)'.*N2(zzB).*Ab0 ...   % stimulated emission term
            + aseTerm(iiS)'.*N2(zzB)./s.Nt ...         % spontaneous emission term
            - alpha*Ab0;                           % attenuation term
        %   centered-space discretized spatial differential
        dAbdz = dzSignB*(Ab1 - Ab0)/(dzSignB*s.dz);
        %   update by forward-marching time differential
        Ab1new = Ab1 + vg*dt*(de - dAbdz);
        Ab1new(Ab1new<0) = 0; % disallow negative values
        Y(zzB,cl.ASEb) = Ab1new;
        %   note change at this point (checking convergence)
        dAbdt = abs(10*log10(Ab1new./Ab1));  
        dAbdt_max(zzB) = max(dAbdt);
        if isnan(dAbdt_max(zzB)) || isinf(dAbdt_max(zzB)), dAbdt_max(zzB) = 100; end           
        
    end
    %%% END POWER PROPAGATION CALCULATIONS FOR THIS UPDATE %%% 
        
    % Calculate inversion based on new fiber conditions
    R12 = Y(:,cl.P)*sig12(iiP)./(h*v(iiP)*Ap/s.Gamma); % Ap already accounts for pump waveguide width; 
                                                        % s.Gamma corrects to account for some light propagating outside bounds of waveguide
    W12 = Y(:,cl.Sf)*sig12(iiS)./(h*v(iiS)*As/s.Gamma);
    W12ASEf = sum( (Y(:,cl.ASEf).*repmat(sig12(iiASE)',[Nz,1])./(h*repmat(v(iiASE)',[Nz 1])*As/s.Gamma)) ,2);
    W12ASEb = sum( (Y(:,cl.ASEb).*repmat(sig12(iiASE)',[Nz,1])./(h*repmat(v(iiASE)',[Nz 1])*As/s.Gamma)) ,2);
    R21 = Y(:,cl.P)*sig21(iiP)./(h*v(iiP)*Ap/s.Gamma);
    W21 = Y(:,cl.Sf)*sig21(iiS)./(h*v(iiS)*As/s.Gamma);
    W21ASEf = sum( (Y(:,cl.ASEf).*repmat(sig21(iiASE)',[Nz,1])./(h*repmat(v(iiASE)',[Nz 1])*As/s.Gamma)) ,2);
    W21ASEb = sum( (Y(:,cl.ASEb).*repmat(sig21(iiASE)',[Nz,1])./(h*repmat(v(iiASE)',[Nz 1])*As/s.Gamma)) ,2);
    A21 = ones(Nz,1)./s.tau;
    N2 = s.Nt * (R12+W12+W12ASEf+W12ASEb) ./ ...
        ((R12+W12+W12ASEf+W12ASEb) + (R21+W21+W21ASEf+W21ASEb+A21));    
    N1 = s.Nt - N2;
    
    if debugFlag
%         [z' 1000*Y(:,cl.Pf) 1000*Y(:,cl.Sf) 1000*sum(Y(:,cl.ASEf),2) 1000*sum(Y(:,cl.ASEb),2) N2/s.Nt N1/s.Nt]
        [z' 1000*Y(:,cl.P) 1000*Y(:,cl.Sf) 1000*sum(Y(:,cl.ASEf),2) 1000*sum(Y(:,cl.ASEb),2) N2/s.Nt N1/s.Nt]
%         plot(z*100,1000*Y(:,cl.Pf),z*100,1000*Y(:,cl.Sf),z*100,1000*sum(Y(:,cl.ASEf),2),z*100,1000*sum(Y(:,cl.ASEb),2));
        plot(z*100,1000*Y(:,cl.P),z*100,1000*Y(:,cl.Sf),z*100,1000*sum(Y(:,cl.ASEf),2),z*100,1000*sum(Y(:,cl.ASEb),2),z*100,1000*s.Pp*N2/s.Nt);
        drawnow;
    end

    % Converge
    if n>1 && max([max(dPdt) max(dSfdt) max(dAfdt_max) + max(dAbdt_max)])<maxChangeAllowed, break; end
    
    %   For each location on fiber, take the last known value and calculate the
    %   differential
    %   Replace all old values with new
    %   Check to see if convergence has occurred

    n = n + 1;                      % count iterations performed
end


% Check analytical solution (neglecting ASE)
if debugFlag
    Ppout = Y(end,cl.Pf);
    Ppin = Y(1,cl.Pf);
    Psout = Y(end,cl.Sf);
    Psin = Y(1,cl.Sf);
    alpha_p = s.GammaP * s.Nt * sig12(iiP);
    Pp_IS = (h*v(iiP)) * (Ap/s.GammaP) * (1/s.tau) * (1/(sig12(iiP)+sig21(iiP)));
    Pp_CS = (h*v(iiP)) * (Ap/s.GammaP) * (1/s.tau) * (1/(sig12(iiP)+sig21(iiS)));
    alpha_s = s.GammaS * s.Nt * sig12(iiS);
    Ps_IS = (h*v(iiS)) * (As/s.GammaS) * (1/s.tau) * (1/(sig12(iiS)+sig21(iiS)));
    Ps_CS = (h*v(iiS)) * (As/s.GammaS) * (1/s.tau) * (1/(sig12(iiP)+sig21(iiS)));
    resP = Ppout - Ppin*exp(-alpha_p*s.L + (Ppin-Ppout)/Pp_IS + (Psin-Psout)/Pp_CS)
    resS = Psout - Psin*exp(-alpha_s*s.L + (Psin-Psout)/Ps_IS + (Ppin-Ppout)/Ps_CS)
    Psout = Psout - resS;
end



if n>maxIterationsAllowed, 
    fprintf('*** Failed to converge after %.2f seconds. ***\n',toc);
    Z = z';
    Yout = zeros(length(Z),5);
    Psout = 0;
    s.displayFlag = false;
    s.displayASEFlag = false;    
else
    fprintf('*** Done after %.2f seconds and %d iterations. ***\n',toc,n);
    % organize output
    Y = Y;                       % Y = [Pf, Pb, Plam1f ... PlamSf ... Plamkf, Plam1b ... PlamSb ... Plamkb], Nz x 2Nlam
    Z = z';                       % n points on z line
    n2 = N2/s.Nt;            % calculate final upper level population
    switch lower(s.direction(1))
        case {'f','r','b'} % forward, reverse, or bidirectionally pumped
            Yout = [Y(:,cl.Pf)+Y(:,cl.Pb) ...       % output column 1: sum of pump powers at all points
                Y(:,cl.Sf) ...              % output column 2: forward signal power
                sum(Y(:,cl.ASEf),2) ...     % output column 3: spectrally integrated forward ASE
                sum(Y(:,cl.ASEb),2)+Y(:,cl.Sb) ... % output column 4: spectrally integrated backward ASE including backward at signal wavelength
                n2];                 % output column 5: inversion ratio (0 to 1)
            Psout = Yout(end,2); % signal output power in W
            PaseF = Yout(end,3);
            PaseB = Yout(1,4);
        case '2' % forward-pumped double pass ... trouble
            Yout = [Y(:,cl.Pf)+Y(:,cl.Pb) ...       % output column 1: sum of pump powers at all points
                Y(:,cl.Sb) ...              % output column 2: backward signal power (this is really the output power)
                sum(Y(:,cl.ASEf),2) ...     % output column 3: spectrally integrated forward ASE
                sum(Y(:,cl.ASEb),2) ...     % output column 4: spectrally integrated backward ASE including backward at signal wavelength
                n2 ...                      % output column 5: inversion ratio (0 to 1)
                Y(:,cl.Sb)];                % output column 6: forward signal power (in case you want it ...)
            Psout = Yout(1,2); % signal output power, W
            PaseF = Yout(end,3);
            PaseB = Yout(1,4);
    end
end

%% Plot

% power plot
if s.plotFlag
    
    hFig = figure('windowstyle','docked');
    Zp = Z*100;                         % [cm] length along fiber
    
    switch(lower(s.mode(1)))
        case 'p'
            % decide if it should be mW or W
            if s.Pp>2,
                ylabStr = 'Pump/Signal Power (W)';
                yr1 = [0 max([s.Pp s.Ps])*1.1];
            else
                ylabStr = 'Pump/Signal Power (mW)';
                Yout(:,1:4)=Yout(:,1:4)*1000;
                yr1 = [0 max([s.Pp s.Ps])*1100];
            end
            % make plot
            [ax h1 h2] = plotyy(Zp,Yout(:,1),Zp,Yout(:,5)); hold on; % pump and inversion ratio
            h3 = plot(ax(1),Zp,Yout(:,2)*10^(-.11),'-or','LineWidth',2); grid on; % signal power
            h4 = plot(ax(1),Zp,Yout(:,3)*10^(-.11),'-ok','LineWidth',2); grid on; % forward ASE
            h5 = plot(ax(1),Zp,Yout(:,4),'-xk','LineWidth',2); grid on; % backward ASE
            %h6 = plot(ax(1),Zp,Y(:,cl.Sb),'-xk','LineWidth',2); grid on; % extra signal power (if desired)
            legVect = [h1 h3 h4 h5 h2]; % pump, signal, ase f, ase b, inversion
            % axis niceties
            set(h1,'LineWidth',2,'Color','b','Marker','o');
            set(h2,'LineWidth',2,'Marker','o');
            set(ax(1),'FontSize',14,'LineWidth',3);
            set(ax(2),'FontSize',14,'LineWidth',3);
            yr2 = [0 min([1.1*max(n2),1])];
            set(ax(1),'YTick',linspace(yr1(1),yr1(2),8)); ylim(yr1);
            set(ax(2),'YTick',linspace(yr2(1),yr2(2),8),'YLim',yr2);
            
            % labels
            xlabel('Fiber Length (cm)');
            ylabel(ylabStr);
            ylabel(ax(2),'Excitation Ratio (N_{2}/N_{tot})');
            set(ax(2),'XTick',[]);
            titleStr = sprintf('%s Power in Fiber -- Forward Pumped',ampName);
            h = title(titleStr); set(h,'FontSize',16);
            % legend
            if s.Pp>2, pumpStr = {sprintf('P_{P}^{+} %dW @ %dnm',ceil(Y(1,1)),round(s.lamP*1e9)),sprintf('P_{P}^{-} %dW @ %dnm',ceil(Y(end,2)),round(s.lamP*1e9))};
            else pumpStr = {sprintf('P_{P}^{+} %dmW @ %dnm',ceil(Y(1,1)*1000),round(s.lamP*1e9)),sprintf('P_{P}^{-} %dmW @ %dnm',ceil(Y(end,2)*1000),round(s.lamP*1e9))};
            end
            legend(legVect, ...
                {sprintf('P_{P} %dmW @ %dnm',round(s.Pp*1000),round(s.lamP*1e9)), ... % pump
                sprintf('P_{S} %dmW @ %dnm',round(Psout*1000),round(s.lamS*1e9)), ... % signal
                sprintf('P_{ASE}^{+} %.1fmW',PaseF*1000), ... % ase forward
                sprintf('P_{ASE}^{-} %.1fmW',PaseB*1000), ... % ase backward
                'Excitation Ratio'},'Location','Best') % excitation ratio
        case 'g'
            % gain values
            Gpf = 10*log10(Y(:,cl.Pf)./s.Pp);     % [dB] pump gain forward
            Gpb = 10*log10(Y(:,cl.Pb)./s.Pp);     % [dB] pump gain backward
            Gsf = 10*log10(Y(:,cl.Sf)./s.Ps);     % [dB] signal gain forward
            Gsb = 10*log10(Y(:,cl.Sb)./s.Ps);     % [dB] signal gain backwards
            Gaf = 10*log10(sum(Y(:,cl.ASEf),2)./s.Ps); % [dB] forward ASE level relative to signal
            Gab = 10*log10(sum(Y(:,cl.ASEb),2)./s.Ps); % [dB] forward ASE level relative to signal
            % plot
            [ax h1 h2] = plotyy(Zp,Gpf,Zp,n2); hold on;
            h3 = plot(ax(1),Zp,Gpb,'-xb','LineWidth',2); grid on;
            h4 = plot(ax(1),Zp,Gsf,'-or','LineWidth',2); grid on;
            h5 = plot(ax(1),Zp,Gsb,'-xr','LineWidth',2); grid on;
            h6 = plot(ax(1),Zp,Gaf,'-ok','LineWidth',2); grid on;
            h7 = plot(ax(1),Zp,Gab,'-xk','LineWidth',2); grid on;
            yr1 = [-10 25];
            ylabStr = 'Pump/Signal Gain (dB)';
            legVect = [h1 h4 h6 h7 h2]; % pump, signal, inversion ... these will normally be all that appear on the legend
            % axis niceties
            set(h1,'LineWidth',2,'Color','b','Marker','o');
            set(h2,'LineWidth',2,'Marker','o');
            set(ax(1),'FontSize',14,'LineWidth',3);
            set(ax(2),'FontSize',14,'LineWidth',3);
            yr2 = [0 min([1.1*max(n2),1])];
            set(ax(1),'YTick',linspace(yr1(1),yr1(2),8)); ylim(yr1);
            set(ax(2),'YTick',linspace(yr2(1),yr2(2),8),'YLim',yr2);
            
            % labels
            xlabel('Fiber Length (cm)');
            ylabel(ylabStr);
            ylabel(ax(2),'Excitation Ratio (N_{2}/N_{tot})');
            set(ax(2),'XTick',[]);
            titleStr = sprintf('%s Single Pass Gain -- Forward Pumped',ampName);
            h = title(titleStr); set(h,'FontSize',16);
            % legend
            if s.Pp>2, pumpStr = {sprintf('P_{P}^{+} %dW @ %dnm',ceil(Y(1,1)),round(s.lamP*1e9)),sprintf('P_{P}^{-} %dW @ %dnm',ceil(Y(end,2)),round(s.lamP*1e9))};
            else pumpStr = {sprintf('P_{P}^{+} %dmW @ %dnm',ceil(Y(1,1)*1000),round(s.lamP*1e9)),sprintf('P_{P}^{-} %dmW @ %dnm',ceil(Y(end,2)*1000),round(s.lamP*1e9))};
            end
            legend(legVect,{sprintf('P_{P}^{+} %dmW @ %dnm',ceil(Y(1,iiP)*1000),round(s.lamP*1e9)), ...
                sprintf('P_{S}^{+} %dmW @ %dnm',ceil(1000*Psout),round(s.lamS*1e9)), ...
                sprintf('P_{ASE}^{+} %dmW',ceil(PaseF*1000)), ...
                sprintf('P_{ASE}^{-} %dmW',ceil(PaseB*1000)), ...
                'Excitation Ratio'},'Location','BestOutside')
        case 'a'
            % get plotting values
            xsurf = Z; % spatial mesh along fiber
            ysurf = lam([2:(iiS-1),(iiS+1):end])'; % ASE wavelengths
            zsurf = Y(:,cl.ASEf) + Y(:,cl.ASEb); % powers in ASE wavelengths
            zsurf = 1000*zsurf;                 % power in mW
            zlabelStr = 'Power (mW)';
            % prepare for plotting
            xsurf = 100*repmat(xsurf,[1,size(zsurf,2)]); % length in cm
            ysurf = 1e9*repmat(ysurf,[size(zsurf,1),1]); % wavelength in nm
            % make surface plot
            surf(xsurf,ysurf,zsurf);
            % labels
            hT = title('Distribution of ASE Power in Amplifier'); set(hT,'fontsize',18);
            set(gca,'fontsize',16,'linewidth',2);
            xlabel('Fiber Length (cm)');
            ylabel('Wavelength (nm)');
            zlabel(zlabelStr);
            % add lines from pump and signal as desired
            shg
            plotMoreChoice = questdlg('Do you want the pump and signal plotted, too?', ...
                'Plot Pump and Signal', ...
                'Both','Just Signal','Neither','Neither');
            switch plotMoreChoice
                case 'Both'
                    hold on;
                    plot3(100*Z,1e9*s.lamP*ones(size(Z)),1000*(Y(:,cl.Pf)+Y(:,cl.Pb)),'-ob'); % pump line
                    plot3(100*Z,1e9*s.lamS*ones(size(Z)),1000*Y(:,cl.Sf),'-or'); % signal line
                    legend('ASE','Pump','Signal');
                    colorbar
                case 'Just Signal'
                    hold on;
                    plot3(100*Z,1e9*s.lamS*ones(size(Z)),1000*Y(:,cl.Sf),'-or'); % signal line
                    legend('ASE','Signal');
                    colorbar
                case 'Neither'
            end
    end
    
else
    hFig = [];
end

%% OUTPUT %% 

if isfield(s,'output')
    if s.output == 'z'
        data{1} = z;
        data{2} = Yout(:,2);
        data{3} = Yout(:,1);
    else
        data = Psout;
    end
else
    data = Psout;
end


end

%%% END MAIN FUNCTION %%%

function aseTerm = GetASETerm(GammaS,sig21,lam,h,c,dlam,Ntot)

aseTerm = 2*dlam*GammaS.*sig21.*(h*(c^2)./(lam.^3))*Ntot;

end