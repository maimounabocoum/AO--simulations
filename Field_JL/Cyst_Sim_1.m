%%%%%%% Simu Field II 02_04_2013 %%%%%%%%%%%%%%%%

%%%% Init step%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
field_init(0);

% height = 5/1000;
% width = 2/1000;
% kerf_x = width/10;%/20;
% kerf_y = height/4;%/8;
% no_ele_x = 10;
% no_ele_y = 4;
% focus = [0 0 40]/1000;
% 
% enabled = ones(no_ele_x,no_ele_y);
% 
% Th = xdc_2d_array(no_ele_x,no_ele_y,width,height, ...
%                     kerf_x,kerf_y,enabled, 1,1,focus);
% xdc_show(Th,'elements');


%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%%
f0 = 4e6;               % Transducer center frequency [Hz]
fs = 200e6;             % Sampling frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
width = 2/1000;         % Width of element
element_height = 5;     % Height of element
kerf = width/20;
Z_focus = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Z_focus];   % Fixed focal point [m]
N_elements = 64;       % Number of elements in the transducer
N_active = 38;          % Active elements in the transducer

% Set the sampling frequency

set_sampling(fs)

% Generate aperture for emission/reception

aperture = xdc_linear_array(N_elements, width, element_height, kerf, 1, 5, focus);

% Set the impulse response and the excitation of the aperture

%%%%%%%% Emission Signal %%%%%%%%%%%%%%%%%%%%%%%

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (aperture, impulse);

% RI = MakeRI_Remote(f0,fs,50);
% Tpulse = length(RI)/fs;
excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(aperture,excitation);

% Load the computer phantom

[phantom_positions, phantom_amplitudes] = cyst_phantom(100);

% Do linear array imaging

no_lines = N_elements-N_active+1;
dx=width;

% Pre-allocate some storage

image_data = zeros(1,no_lines);

for i=1:no_lines
    i
    % Find position for imaging
    x=(i-1-no_lines/2)*dx;
    
    % Set the focus for this direction
    
    xdc_center_focus(aperture,[x 0 0]);
    xdc_focus(aperture,0,[x 0 Z_focus]);
    xdc_center_focus(aperture,[x 0 0]);
    xdc_focus(aperture,0,[x 0 Z_focus]);

    % Set the active elements using the apodization
    
    apo = [zeros(1,i-1) hamming(N_active)' zeros(1,N_elements-N_active-i+1)];
    xdc_apodization(aperture,0,apo);
    xdc_apodization(aperture,0,apo);
    
    % Calculate the received response
    
    [v, t1] = calc_scat_multi(aperture,aperture,phantom_positions, phantom_amplitudes);
    
    % Store the result
    
    image_data(1:max(size(v)),i)=v;
    times(i)=t1;
end

% Free space for apertures
xdc_free(aperture)

% image

min_sample = min(times)*fs;
for i=1:1:no_lines
    rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
    env(1:size(rf_env,1),i) = rf_env;
end

env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
figure;
image(x*1000,depth*1000,env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
coloramp(gray(128))
title('Image of cyst phantom (60dB dynamic range)')









