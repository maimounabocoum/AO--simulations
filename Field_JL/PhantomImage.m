clear all
close all
clc

field_init(0);
CurrentDirectory = cd;

% Set initial parameters for Scatterers
f0=5.5e6;           % Transducer center frequency [Hz]
fs=200e6;           % Sampling frequency [Hz]
c=1540;             % Speed of sound [m/s]
lambda=c/f0;        % Wavelength [m]

N_elements = 1;
focus = [0 0 40]/1000;
R_elem = 3.5/1000;
el_size = 0.4/1000;

Th = xdc_concave(R_elem, focus(3), el_size);

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, impulse);
figure; plot(excitation)

% Scatterers Definition

x_phantom_size = 3/1000;
y_phantom_size = 3/1000;
z_phantom_size = 5/1000 - 1/1000;
z_phantom_start = 38/1000;
N_scat = 2000;

x = (rand(N_scat,1)-0.5)*x_phantom_size;
y = (rand(N_scat,1)-0.5)*y_phantom_size;
z = rand(N_scat,1)*z_phantom_size + z_phantom_start;

phantom_amplitudes = randn(N_scat,1);
phantom_positions = [x y z];
clear x y z

save(['phantom_amplitudes'],'phantom_amplitudes');
save(['phantom_positions'],'phantom_positions');
disp('Data saved')


