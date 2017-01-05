%% Example Title
% Summary of example objective

%% Example Title
clearvars
close all
addpath('..\Field_II')
field_init(0);

%%  creation of a single transducer :
% Set initial parameters

f0=3e6; % Transducer center frequency [Hz]
f1=1e6; % Test frequency 1 [Hz]
f2=10e6; % Test frequency 2 [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength [m]
height=10/1000; % Height of element [m]
width=lambda; % Width of element [m]
kerf=width/4; % Distance between transducer elements [m]
N_elements=2; % Number of elements
no_sub_x=1; % Number of sub-divisions in x-direction of elements.
no_sub_y=1; % Number of sub-divisions in y-direction of elements.
focus=[0 0 80]/1000; % Initial electronic focus

% initialization of sound velocity:
%set_field ('c',c); % in m/s

% Set the sampling frequency
%set_sampling(fs);

% Define the transducer as a pointer
Th = xdc_linear_array (N_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
%show_xdc (Th) 


% Set the waveforms
waveform1=sin(2*pi*f1*(0:1/fs:4/f1));
element_no=1;
ele_waveform (Th, element_no, waveform1);
waveform2=sin(2*pi*f2*(0:1/fs:4/f2));
element_no=2;
ele_waveform (Th, element_no, waveform2);


% difnition of the line over which the impulse response in calculated:

[h,t] = calc_h (Th,[0 0 60]/1000);

figure;
plot((0:length(h)-1)/fs+t,h,'color','red','marker','o')
ylabel('Response in m/s')
xlabel('Time \mu s')
 
length(h)
field_info;

field_end;



