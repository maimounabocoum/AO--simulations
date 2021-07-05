%%%%% initialize programm FIELD II %%%%%%%%%%%%
% Maïmouna BOCOUM - last edited versions 03-04-2019
% note : all parameters are defined in SI units

%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%

% set_field('show_time',1);

param.f0 = 3e6;                     % Transducer center frequency [Hz]
param.fs = 100e6;                   % Sampling frequency in FIELDII[Hz]
param.fs_aq  = 10e6;                % Sampling frequency of the photodiode [Hz]
param.Noc = 4 ;                     % Number of optical cycles
param.c = 1540;                     % Speed of sound [m/s]
param.lambda = param.c/param.f0;    % Wavelength [m]
param.element_height= 6/1000;       % Height of element [m] 6
param.width = 0.2/1000;             % Width of element [m] - 0.11 for 15MhZ probe
param.kerf = 0/1000;                % Distance between transducer elements [m]
param.N_elements = 192;             % 192; % Number of elements for SL10-2 probe
param.X0 = -50/1000  ;              % position min of effective probe shooting (center probe = 0mm)
param.X1 =  50/1000 ;               % position max of effective probe shooting (center probe = 0mm)
param.Rfocus = 35/1000;             % Static Elevation focus
param.attenuation = 0;              % en db/cm/MHz (not yet implemented)- in water : 0.3 ?
param.no_sub_x = 1;
param.no_sub_y = 11; % for designed probes, you should put a value > 2 for proper calculation (10 is good!)


param.farfield  = param.width^2/(4*param.lambda); 

param.detection = {'camera','photodiode'};        % 'photodiode': return flux over time / 'camera': returned fluence over camera integration time
param.tau_c     = 20e-6;        % camera intergration type for holography detection (advise: leave it at 20us )
param.Trigdelay = 30e-6;        % accounts for time fo flight delay before camera triggers

% tips: you can evaluate tau_c and Trigdelay using the "visualize the
% field" bloc adjusting parameters in:
% CurrentExperiement.ShowFieldCorrelation()

%% type of focalization to apply for the virtual experiment :
% OF : 'Focused Waves'
% OP : 'Plane Waves'
% JM : 'Jean-Michel continuoous waves' (not implemented yet)
% OS : 'Structured Waves ' 

param.FOC_type    = 'JM'; 
param.Bascule     = 'off';             % parameter for JM with / without Talbot Effect
param.focus       = 23/1000;           % Initial electronic focus     - only active in OF mode
param.angles      = 0*pi/180;          % Line Vector Angular scan     - only active in OP and OS mode 


% k0 = (1/1e-3) is the smapling frequence for the decimation
% k0 = (1/(param.N_elements*param.width)) is the smapling frequence for the decimation

param.df0x = (1/(param.N_elements*param.width));  % 24.39; - only active in OP and OS mode 
param.decimation  = 10;  % decimation list of active actuators   - only active in OS mode 
% decimation definition : 
% activeElements are indexed by 
% mod( (1:N_elements) - ElmtBorns(1) , 2*decimation ) ;

param.NbZ         = 8;                          % 8; % Nb de composantes de Fourier en Z, 'JM'
param.NbX         = 3;                          % 20 Nb de composantes de Fourier en X, 'JM'
param.phase       = 0;                          % phases i 2pi unit for 'JM'
param.nuZ0 = 1/( (param.c)*20e-6 );             % Pas fréquence spatiale en Z (en mm-1)
param.nuX0 = 1/(param.N_elements*param.width);  % Pas fréquence spatiale en X (en mm-1) 



param.Activated_FieldII = 1 ;     % 0 to generate field by yourself - 1 FIELDII simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Simulation BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation box initialization : 

    param.Xrange = [-15 15]/1000;     % [-15 15]/1000;     % in m [-15 15]
    param.Yrange = 0/1000;            % [-0.1 0.1]/1000 ; (not implemented yet)
    param.Zrange = [20 50]/1000;      % [0.5 40]/1000;       % simulation JM : [5 40]/1000;

    param.Nx = 200;             % number of interpolating points along Xrange
    param.Ny = 1;               % number of interpolating points along Yrange
    param.patternRep = 4;       % number of times the 40us main pattern is repeted (minimum = 1) 
    % in order to match fs_aq(Hz) along Zrange , and 
    % unshures Nz >=1
    param.Nz = 160;%max( 1 , ceil ( param.fs_aq * (abs(param.Zrange(2) - param.Zrange(1)))/(param.c) ) ); % do not edit
%% definition of laser beam
    
% waist of diffuse IR laser beam
param.w0 = [10 10]/1000 ;             % specify the center of the gaussian beam.
param.center = [0 0 32.5]/1000 ;      % specify the center of the gaussian beam.    
             
                                    % if this value is commented, 
                                    % the beam is by defaukt center on the
                                    % simulation box
%% absorbers positions :
    % fringes : modulation of intensity in direction given by Position
    param.phantom.Positions = [0 0 31.5 ; 0 0 33.5]/1000;  % [x1 y1 z1; x2 y2 z2 ; ....] aborbant position list
    param.phantom.Sizes     = [1 ; 1]/1000;             % dimension in all direction [dim ; dim ; ...]    param.phantom.Types = {'gaussian','gaussian'} ;         % available types exemple : { 'square', 'gaussian', ...}
    param.phantom.Types = {'gaussian','gaussian'};
    
%% Probe defintion :
% Set the sampling frequency
if param.Activated_FieldII == 1 
field_debug(0);
set_sampling(param.fs);
set_field('c',param.c);
set_field('Freq_att',param.attenuation*100/1e6);
set_field('att',0*2.6*100);
set_field ('att_f0',0*param.f0); 
set_field('use_att',1);
end

% screen field parameters :

