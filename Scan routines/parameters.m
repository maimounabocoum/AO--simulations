%%%%% initialize programm FIELD II %%%%%%%%%%%%

%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%

% set_field('show_time',1);

param.f0 = 6e6;                     % Transducer center frequency [Hz]
param.fs = 100e6;                   % Sampling frequency using in FIELDII[Hz]
param.fs_aq  = 10e6;                % Sampling frequency for the experiement [Hz]
param.Noc = 2 ;                     % Number of optical cycles
param.c = 1540;                     % Speed of sound [m/s]
param.lambda = param.c/param.f0;    % Wavelength [m]
param.element_height= 6/1000;       % Height of element [m] 6
param.width = 0.2/1000;             % Width of element [m] - 0.11 for 15MhZ probe
param.kerf = 0/1000;                % Distance between transducer elements [m]
param.N_elements = 192;            % 192; % Number of elements
param.X0 = -19.2/1000  ;               % prosition min of effective probe shooting
param.X1 =  19.2/1000 ;                % prosition min of effective probe shooting
param.Rfocus = 35/1000;             % Static Elevation focus
param.attenuation = 0;              % en db/cm/Mhz
param.no_sub_x = 1;
param.no_sub_y = 10; %for designed probes, you should put a value > 2 for proper calculation (10 is good!)
param.farfield = param.width^2/(4*param.lambda); 

%% type of focalization to apply for the virtual experiment :
% OF : 'Focused Waves'
% OP : 'Plane Waves'
% JM : 'Jean-Michel continuoous waves'
% OS : 'Structured Waves ' 
param.FOC_type = 'OS'; 

param.focus       = 10/1000;              % Initial electronic focus [m,m,m]      - only active in OF mode
param.angles      = (-90:90)*pi/180;    % Angular scan [m,m,m]                  - only active in OP and OS mode 
% k0 = (1/1e-3) is the smapling frequence for the decimation
% k0 = (1/(param.N_elements*param.width)) is the smapling frequence for the decimation

param.df0x = (1/(param.N_elements*param.width)) ;
param.decimation  = [10];  % decimation list of active actuators   - only active in OS mode 
% decimation definition : 
% activeElements are indexed by 
% mod( (1:N_elements) - ElmtBorns(1) , 2*decimation ) ;

param.NbZ         = 10;                  % 8; % Nb de composantes de Fourier en Z, 'JM'
param.NbX         = -10:10;               % 20 Nb de composantes de Fourier en X, 'JM'
 



param.Activated_FieldII = 1 ;     % 0 to generate field by yourself

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Simulation BOX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation box initialization : 

    param.Xrange = [-15 15]/1000;  % in m
    param.Yrange = 0/1000;        % [-0.1 0.1]/1000;
    param.Zrange = [2 35]/1000;   % in m

    param.Nx = 150;
    param.Ny = 1;
    % in order to match the number of point in Z direction , and 
    % unshures Nz >=1
    param.Nz = max( 1 , ceil ( param.fs_aq * (abs(param.Zrange(2) - param.Zrange(1)))/(param.c) ) ); % do not edit
% waist of diffuse IR laser beam
param.w0 = [9 9]/1000 ; % in m 
param.center = [0 0 18.5]/1000 ;      
             % specify the center of the gaussian beam.
                                    % if this value is commented, 
                                    % the beam is by defaukt center on the
                                    % simulation box
    %% abosbers positions :
    % fringes : modulation of intensity in direction given by Position
    param.phantom.Positions = [-1.5*0 0 18.5 ; 1000 0 18.5]/1000; % [x1 y1 z1; x2 y2 z2 ; ....] aborbant position list
    param.phantom.Sizes     = [0.8 ; 0.8]/1000;          % dimension in all direction [dim ; dim ; ...]
    param.phantom.Types = {'gaussian','gaussian'} ;   % available types exemple : { 'square', 'gaussian', ...}
    
    
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

