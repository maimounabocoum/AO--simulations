%%%%% initialize programm FIELD II %%%%%%%%%%%%

%%%%%%%%% set parameters %%%%%%%%%%%%%%%%%%%%
field_debug(0);
% set_field('show_time',1);

f0=6e6; % Transducer center frequency [Hz]
fs=200e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
element_height= 6/1000; % Height of element [m]
width=0.2/1000; % Width of element [m]
kerf= 0; % Distance between transducer elements [m]
N_elements = 128;%128; % Number of elements
ActiveList = 1:128; % index of active elements
focus = [0 0 40]/1000; % Initial electronic focus
Rfocus = 70/1000; % Elevation focus
attenuation = 0;         % en db/cm/Mhz
no_sub_x = 5;
no_sub_y = 10;
farfield = width^2/(4*lambda); 


%% Probe defintion :
% Set the sampling frequency
set_sampling(fs);
set_field('c',c);
set_field('Freq_att',attenuation*100/1e6);
set_field('att',0*2.6*100);
set_field ('att_f0',f0); 
set_field('use_att',1);