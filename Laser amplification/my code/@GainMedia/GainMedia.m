classdef GainMedia
    %GAINMEDIA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % rod caracteristics
        Rod
        L  = 10e-2;     % crystal length in m
        w0  = 100e-6;   % active surface 
        
        tau % lasing fluorescent lifetime
        gamma 
        N0  % laser density
    end
    
    methods
        
        
        function obj = GainMedia(Rod)
            %GAINMEDIA Construct an instance of this class
            %   Detailed explanation goes here
            obj.Rod = Rod ;
switch Rod
    case 'Ruby'
        obj.tau   = 3e-3;         % fluorescent time
        obj.lambda_e = 694.4e-9;
        nu_e = c/lambda_e ;
        sigma_e = 475e-9*1e-4; % emission cross section m2 at 750nm
        obj.Is = (Ee)/(tau*sigma_e); % Saturation intensity W/m2
        % pumpinp with flashlamp
        obj.eta = 0.4; % absorption slope
        obj.lambda_p = 475e-9; % pump wavelength        
        nu_p = c/lambda_p;
        obj.Ep = h*nu_p ;
        obj.N0 = 1.58e19*1e6 ; % Cr3+ concentration in cm^{-3}
        
    case 'Alexandrite'
        tau   = 260e-6;         % fluorescent time at 298K
        lambda_e = 750e-9;
        nu_e = c/lambda_e ;
        Ee = h*nu_e ;           % photon energy
        sigma_e = 1e-20*1e-4; % emission cross section m2 at 750nm
        GainBW = 50e-9; % gain Bandwith m
        %https://books.google.fr/books?id=_F4hAwAAQBAJ&pg=PA564&lpg=PA564&dq=small+signal+gain+in+alexandrite&source=bl&ots=PO4yBTTSNg&sig=ACfU3U2xkBBNMXwOiQdgyRImEOMXa6XMXw&hl=fr&sa=X&ved=2ahUKEwier_zs2prpAhUHLBoKHcswDT4Q6AEwEnoECAgQAQ#v=onepage&q=small%20signal%20gain%20in%20alexandrite&f=false
        Is = (Ee)/(tau*sigma_e); % Saturation intensity W/m2
                % pumpinp with flashlamp
        eta = 0.4; % absorption slope
        lambda_p = 630e-9; % pump wavelength        
        nu_p = c/lambda_p;
        Ep = h*nu_p ;
        % small gain coefficient 4-20 /m
        % inversion density 6e24 /m3
        N0 = 3.52e19*1e6 ; % 0.1%-doping concentration in cm^{-3}
    case 'Nd:YVO4'
        obj.tau   = 90e-6;             % fluorescent time
        obj.gamma = 1;                 % degenerency ratio
        obj.lambda_e = 1064e-9;
        obj.nu_e = c/lambda_e ;
        obj.Ee = h*nu_e ;           % photon energy
        obj.sigma_e = 25e-19*1e-4;  % emission cross section m2 at 1064nm
                                % https://www.unitedcrystals.com
        obj.GainBW = 1e-9; % gain Bandwith m
        obj.Es = 1e-4*(h*nu_e)/(gamma*sigma_e); % saturation fluence J/cm2
        obj.Is = (Ee)/(tau*sigma_e); % Saturation intensity W/m2
        % pumpinp at 946nm
        obj.N0 = 1.25e20*1e6 ;          % doping concentration cm^{-3}
                                    % 1%-doping concentration in cm^{-3}
                                    % https://www.rp-photonics.com/doping_concentration.html
        obj.eta = 0.48;                 % absorption slope
                                    % http://www.pmoptics.com/neodymium_doped_yvo4.html
        obj.sigma_a = 0.5*( 8e-20 + 24e-20 )*1e-4; % absorption cross section m2 at 808nm p.86
        % 30 cm-1 : sigma_a = 24e-20*1e-4 along C-axis
        % 10 cm-1 : sigma_a = 8e-20*1e-4 along A-axis
        obj.AbsBW = 4e-9;               % absortption Bandwith m
        obj.lambda_p = 808e-9;          % pump wavelength        
        obj.nu_p = c/lambda_p;
        obj. Ep = h*nu_p ;
    case 'Nd:YAG'
        obj.tau   = 230e-6;         % fluorescent time
        lambda_e = 1064e-9;
        nu_e = c/lambda_e ;
        Ee = h*nu_e ;           % photon energy
        sigma_e = 2.8e-19*1e-4; % emission cross section m2 at 1064nm
        % https://www.rp-photonics.com/yag_lasers.html
        GainBW = 0.6e-9; % gain Bandwith m
        Es = 1e-4*(h*nu_e)/(gamma*sigma_e); % saturation fluence J/cm2
        Is = (Ee)/(tau*sigma_e); % Saturation intensity W/m2
        
        % pumpinp at 808nm
        N0 = 1.38e20*1e6 ;      % 1%-doping concentration in cm^{-3}
        sigma_a = 4.1e-20*1e-4; % absorption cross section m2 at 808nm p.86
        % eta = 0.32; % absorption slope
        lambda_p = 808e-9; % pump wavelength        
        nu_p = c/lambda_p;
        Ep = h*nu_p ;
        
    case 'TiSa'
        % emission properties
        tau = 3.2e-6;        
        gamma = 1;          % degenerency ratio
        lambda_e = 790e-9;  % pump wavelength
        sigma_e = 4.1e-19*1e-4; % emission cross section m2 at 790nm
        % https://www.rp-photonics.com/titanium_sapphire_lasers.html?s=ak
        nu_e = c/lambda_e;
        Ee = h*nu_e ;       % photon energy
        Es = 1e-4*(h*nu_e)/(gamma*sigma_e) ; % saturation fluence J/cm2
        Is = (Ee)/(tau*sigma_e);             % saturation intensity W/m2
        GainBW = 230e-9;    % gain Bandwith m (p.93)
        
        % pumpinp at 532nm propreties
        eta = 0.34; % absorption slope
        lambda_p = 532e-9; % pump wavelength        
        nu_p = c/lambda_p;
        Ep = h*nu_p ;
        N0 = 4.56e19*1e6 ; % Cr3+ concentration in cm^{-3} -> m^{-3}

       

end


        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

