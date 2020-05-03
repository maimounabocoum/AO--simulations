%% constants
h = 6.62e-34;% J.s
c = 3e8;
%% rod parameters
L  = 1e-2;     % crystal length in m
w0  = 100e-6; % active surface 

%% pump caracteristics
Pp = 10 ; % pump power in Watt
lambda_p = 532e-9; % pump wavelength
nu_p = c/lambda_p;
Ep = h*nu_p ;
n = 1.76 ; % index of refraction
eta = 0.38 ; % slope efficiency
Is = 2.85e9; % saturation intensity W/m2 (Ozawa PRL 2010)
kappa = 8.84e-29; % m2.s
Rp = Pp*eta/(L*w0^2*pi*Ep);

%% steady state

tau_f = 3.2e-6;        
gamma = 1;      % degenerency ratio
lambda_e = 790e-9;
nu_e = c/lambda_e ;
sigma_e = 4.1e-19*1e-4; % emission cross section m2 at 790nm
% https://www.rp-photonics.com/titanium_sapphire_lasers.html?s=ak
Es = 1e-4*(h*nu_e)/(gamma*sigma_e); % saturation fluence J/cm2
GainBW = 230e-9; % gain Bandwith m (p.93)
% pumpin at 532nm


% Pin / Pout vs Pump power in steady state
Pin = 0.1;
Pout = 1: 100;
[PIN, POUT] = meshgrid(Pin,Pout);

PPUMP = (Ep/(kappa*eta))*(pi*w0^2*log(POUT./PIN) + (POUT-PIN)/Is);

figure(1)
plot(PPUMP,Pout)
xlabel('Pump Power')
ylabel('output Power')
%cb = colorbar ;
%ylabel(cb,'W')








