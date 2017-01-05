%%%%%%% Simu Field II 02_04_2013 %%%%%%%%%%%%%%%%

%%%% Init step%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
field_init(0);

% Set initial parameters
f0=4e6;             % Transducer center frequency [Hz]
fs=200e6;           % Sampling frequency [Hz]
c=1540;             % Speed of sound [m/s]
lambda=c/f0;        % Wavelength [m]

N_elements = 1;
width=4/1000;       % Width of element [m]
height=2/1000;      % Height of element [m]
kerf=width/4;       % Distance between tr
Rconvex = 5/1000;
Rfocus = 1/1000;
no_sub_x = 10;      % Warnin: Field2 error, x and y inverted
no_sub_y = 10;      % Warnin: Field2 error, x and y inverted
focus = [0 0 40]/1000;

%Th = xdc_linear_array(N_elements, width, height, kerf, 2, 3, focus);
Th = xdc_convex_focused_array(N_elements, width, height, kerf, Rconvex, Rfocus, no_sub_x, no_sub_y, focus);
xdc_show(Th,'elements');
%dfgdfg
impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, impulse);
figure; plot(excitation)

[phantom_positions, phantom_amplitudes] = phantom(2000);

%%%% ROT et trans
Rotate = [cos(atan(-0.1388)) 0 sin(atan(-0.1388)) ; 0 1 0 ; -sin(atan(-0.1388)) 0 cos(atan(-0.1388))];
Translate = [-5.5/1000 0 0.3831/1000];
for i=1:1:size(phantom_positions,1)
    phantom_positions2(i,:) = phantom_positions(i,:)*Rotate + Translate;
end

% [v,t]= calc_scat_multi(Th,Th,[0 0 20]/1000,1);
%[v,t]= calc_scat_multi(Th,Th,phantom_positions,phantom_amplitudes);
[v,t]= calc_scat_multi(Th,Th,phantom_positions2,phantom_amplitudes);

% min_sample = t*fs;
% 
% rf_env=abs(hilbert([zeros(round(t*fs-min_sample),1); v(:)]));
% env(1:size(rf_env,1)) = rf_env;


%subplot(211)
figure
[N,M]=size(v);
v=v/max(max(v));
for i=1:N_elements
    plot((0:N-1)/fs+t,v(:,i)+i), hold on
end
hold off
title('Individual traces')
xlabel('Time [s]')
ylabel('Normalized response')
%subplot(212)
figure
if size(v,2)>1
    plot((0:N-1)/fs+t,sum(v'))
else
    plot((0:N-1)/fs+t,v(:,i)+i)
end
title('Summed response')
xlabel('Time [s]')
ylabel('Normalized response')




