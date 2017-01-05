%%%%%%% Simu Field II 29_03_2013 %%%%%%%%%%%%%%%%

%%%% Init step%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
field_init(0);

%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%%
R = 10/1000;            % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 1/1000;      % Size of math elements
f0 = 4e6;               % Transducer center frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
No_elements = 1;        % # of elements
times = 0;
P0 = 40e5;              % Pression mesure au foyer en manip

InfoBeam.Xfocus = focus(1);
InfoBeam.Yfocus = focus(2);
InfoBeam.Zfocus = focus(3);

%%%% Field Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
Rho = 1000;             % Density
fs = 200e6;             % Sampling frequency [Hz]
attenuation = 0.6;      % in db/cm/Mhz

%%%%%%%% Transducer Definition %%%%%%%%%%%%%%%%%

Th = xdc_concave(R,Rfocal,ele_size); % Define concave source
xdc_show(Th,'elements');

%%%%%%%% Elements Definition %%%%%%%%%%%%%%%%%%%



%%%%%%%% Emission Signal %%%%%%%%%%%%%%%%%%%%%%%

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
RI = MakeRI_Remote(f0,fs,50);
Tpulse = length(RI)/fs;
xdc_excitation(Th,RI);

%%%%%%%% Apodisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%apodisation = ones(No_elements,1);
%xdc_apodization(Th,times,apodisation');

[scat, star_time] = calc_scat_multi(Th,Th,points,amplitudes); % Procedure for calculating the received signal from a collection of scatterers for each of the elements of the aperture Th





