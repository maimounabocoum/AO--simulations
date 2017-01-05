clear all
close all
clc
field_init(0);

% Set initial parameters
f0=3e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
height=5/1000; % Height of element [m]
width=1/1000; % Width of element [m]
kerf=width/4; % Distance between tr
N_elements = 32;
focus = [0 0 40]/1000;

Th = xdc_linear_array(N_elements, width, height, kerf, 2, 3, focus);

impulse_response = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse(Th,impulse_response);

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation(Th, excitation);

[v,t]=calc_scat_multi(Th,Th,[0 0 20]/1000,1);


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
axis([t t+N/fs 0 M+1])
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
axis([t t+N/fs 0 M+1])