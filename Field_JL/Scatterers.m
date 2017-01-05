%%%%%%% Simu Field II 02_04_2013 %%%%%%%%%%%%%%%%

%%%% Init step%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
field_init(0);
CurrentDirectory = cd;

% Set initial parameters
f0=5.5e6;             % Transducer center frequency [Hz]
fs=200e6;           % Sampling frequency [Hz]
c=1540;             % Speed of sound [m/s]
lambda=c/f0;        % Wavelength [m]

N_elements = 1;
focus = [0 0 40]/1000;
R_elem = 3.5/1000;
el_size = 0.4/1000;

Th = xdc_concave(R_elem, focus(3), el_size);
%Th = xdc_linear_array(N_elements, width, height, kerf, 2, 3, focus);
%Th = xdc_convex_focused_array(N_elements, width, height, kerf, Rconvex, Rfocus, no_sub_x, no_sub_y, focus);

%xdc_show(Th,'elements');

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, impulse);
figure; plot(excitation)

% Scatterers
[phantom_positions, phantom_amplitudes] = phantom(2000);


[name,directory]=uiputfile('phantom_');
if isequal(name,0) || isequal(directory,0);
else
    cd(directory);
    save([name 'positions.mat'],'phantom_positions');
    save([name 'amplitudes.mat'],'phantom_amplitudes');
    disp('Data saved')  
    cd(CurrentDirectory);     
end



