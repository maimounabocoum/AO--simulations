%% Example Title
% Summary of example objective

%% Example Title
clearvars
%%  creation of a single transducer :
% Set initial parameters
f1 = 1e6; % waform frequence in Hz
fs = 100e6; % sampling frequency

c = 1500;
height=5/1000; % Height of element [m]
width=2/1000; % Width of element [m]
kerf=width/4; % Distance between transducer elements [m]
N_elements=1; % Number of elements
focus=[0 0 35]/1000; % Initial electronic focus

% initialization of sound velocity:
set_field ('c',c); % in m/s

% Define the transducer as a pointer
Th = xdc_linear_array (N_elements, width, height, kerf, 2, 3, focus);

% set a waveform :
dt = 1/fs;
tmax = 4e-6; % end of waveform

% figure(1)
 waveform1=sin(2*pi*f1*(0:dt:tmax));
% plot(2*pi*f1*(0:dt:tmax)*1e6,waveform1,'marker','o');
% xlabel('Time in \mu s')
% application of the waveform on the transducer :
ele_waveform (Th,1, waveform1);


% difnition of the line over which the impulse response in calculated:

z = [1:0.1:10]/1000;

VelocityResponse = 0*z; % initialize in m/s
tsample = [0:0.01:10]; % time tabulation in \mu s
[Z, Tsample] = meshgrid(z,tsample);

Line = zeros(length(z),3);
Line(:,:,3) = z;
[h,t] = calc_h (Th,Line);

%HResponse(loop,:) = interp1(((0:length(h)-1)*dt+t)*1e6,h(:,1),tsample,'linear',0);

% hold on
% plot(tsample,HResponse,'color','red','marker','o')
% %plot(((0:length(h)-1)*dt+t)*1e6,h(:,1),'color','red','marker','o')
% ylabel('Response in m/s')
% xlabel('Time \mu s')
% pause(0.5)


clf; hf = figure(1);
for loop = 1:10:length(tsample)
%imagesc(tsample,z*1e3,HResponse)
plot(z*1e3,HResponse(:,loop))
ylim([0 3000])
title(['time = ',num2str(tsample(loop))])
pause(0.1)
end
