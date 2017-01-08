%%% this is a test routine %%%
%% maimouna bocoum 04-01-2017
clear all
clearvars ;

addpath('..\Field_II')
field_init(0);
parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate aperture for emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = ActuatorProbe(N_elements,element_height,width,no_sub_x,no_sub_y,kerf,ActiveList);
Probe = xdc_rectangles(P.rect,[0 0 0], focus);
%Probe = xdc_linear_array (N_elements, width, element_height, kerf,no_sub_x,no_sub_y, focus);
%Probe = xdc_focused_array(N_elements,width,element_height,kerf,Rfocus,no_sub_x,no_sub_y,focus);
%data = xdc_get(Probe,'rect');
%show_xdc (Probe)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting the impulse response field (Green function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emission signal : I think by befault, the green function is calculated
% for a different frequency
t_impulseResponse = (0:1/fs:2/f0);
impulse = sin(2*pi*f0*t_impulseResponse);
impulse=impulse.*hanning(length(impulse))'; 


xdc_impulse (Probe, impulse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of excitation function for actuator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = N_elements;
PositionActuators = [1:N]';
% user defined excitation :
% N number of actuator 
% PositionActuators : Position of active actuator. 

%excitation= Actuators.Excitation(1,:);%sin(2*pi*f0*(0:1/fs:2/f0));
 Noc = 2;
t_excitation = (0:1/fs:Noc*1.5/f0);
excitation =  sin(2*pi*f0*t_excitation);
excitation = excitation.*hanning(length(excitation))';

% delay law :
delays = DelayLaw(width,N_elements,c,40/1000);

%xdc_focus(Probe,0,[0 0 30]/1000);
xdc_focus_times(Probe,-1, delays);
xdc_excitation (Probe, excitation);
%xdc_show(Probe,'elements')


% apodisation :
% apodisation = ones(N_elements,1);
% xdc_apodization(Probe,0,apodisation')
%xdc_focus (Probe,4e-6,[0 0 3]/1000);

%% Simulation box initialization : 
Nx = 50;
Ny = 1;
Nz = 40;

Xrange = [-5 5]/1000; % in m
Yrange = 0;%[-0.1 0.1]/1000;
Zrange = [25 45]/1000; % in m

SimulationBox = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz);

% calculation of the emitted Field :
[h,t] = calc_hp(Probe,SimulationBox.Points());
h = h/max(h(:));
time = t + [0:(size(h,1)-1)]/fs;
% h : each column corresponds to the field calculated for the
% SimulationBox.Points() list

% Screen maximum field for each position:
Field_resized = reshape(h',[Ny,Nx,Nz,length(time)]);
% first dimension now corresponds to the y values.


% initialize frames for movie capture :

%Screen Raw data Result
Hf2 = figure(2);
set(Hf2,'name','raw data returned By calc_hp function')
imagesc(1:Nx*Ny*Nz,time,abs(h./repmat(max(h),size(h,1),1)).^2)
title('image in log scale')
xlabel('Point index')
ylabel('time (\mu s)')
colorbar




xdc_free(Probe);
field_end;
