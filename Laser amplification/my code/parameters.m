%% constants
h = 6.62e-34;% J.s
c = 3e8;
%% parameters sheet
Rod = 'Nd:YAG';
L  = 1e-2;     % crystal length in m
w0  = 100e-6; % active surface 

%% pump caracteristics
Pp = 10 ; % pump power in Watt
lambda_p = 532e-9; % pump wavelength
nu_p = c/lambda_p;
Rp = Pp*eta/(L*w0^2*pi*h*nu_p);

%%

switch Rod
    case 'Nd:YAG'
        tau_f = 230e-6; % fluorescent time
        gamma = 1;      % degenerency ratio
        lambda_e = 1064e-9;
        nu_e = c/lambda_e ;
        sigma_e = 2.8e-19*1e-4; % emission cross section m2 at 1064nm
        % https://www.rp-photonics.com/yag_lasers.html
        GainBW = 0.6e-9; % gain Bandwith m
        Es = 1e-4*(h*nu_e)/(gamma*sigma_e); % saturation fluence J/cm2
        % pumpin at 946nm
        Is = 2.9e7 ; % Saturation intensity W/m2
        N0 = 1.5e20; % cm-3 : total dopant concentration
    
    case 'TiSa'
        tau_f = 3.2e-6;        
        gamma = 1;      % degenerency ratio
        lambda_e = 790e-9;
        nu_e = c/lambda_e ;
        sigma_e = 4.1e-19*1e-4; % emission cross section m2 at 790nm
        % https://www.rp-photonics.com/titanium_sapphire_lasers.html?s=ak
        Es = 1e-4*(h*nu_e)/(gamma*sigma_e); % saturation fluence J/cm2
        GainBW = 230e-9; % gain Bandwith m (p.93)
        % pumpin at 532nm
end










