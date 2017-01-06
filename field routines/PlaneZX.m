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

Probe = xdc_linear_array (N_elements, width, element_height, kerf,no_sub_x,no_sub_y, focus);
%Probe = xdc_focused_array(N_elements,width,element_height,kerf,Rfocus,no_sub_x,no_sub_y,focus);
%show_xdc (Probe)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting the impulse response field (Green function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emission signal : I think by befault, the green function is calculated
% for a different frequency
t_impulseResponse = (0:1/fs:2/f0);
impulse = sin(2*pi*f0*t_impulseResponse);
impulse=impulse.*hanning(length(impulse))'; 

% figure;
% plot(1e6*t_impulseResponse,impulse);
% xlabel('time(\mu s)')
% ylabel('a. u')
% title('Probe impulse Response setting')

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
delays = DelayLaw(width,N_elements,c,10/1000);
figure
plot(delays)

%xdc_focus(Probe,0,[0 0 30]/1000);
xdc_focus_times(Probe,-1, delays);

% figure;
%  plot(1e6*t_excitation,excitation);
%  xlabel('time(\mu s)')
%  ylabel('a. u')
%  title('Probe uniform excitation')


xdc_excitation (Probe, excitation);
%xdc_show(Probe,'focus')
  
  
% RI = MakeRI_Remote(f0,fs,50);
% Tpulse = length(RI)/fs;
% xdc_excitation (Probe, RI);

% apodisation :
% apodisation = ones(N_elements,1);
% xdc_apodization(Probe,0,apodisation')
%xdc_focus (Probe,4e-6,[0 0 3]/1000);
%ele_waveform(Probe,PositionActuators,Actuators.Excitation);

%% Simulation box initialization : 
Nx = 50;
Ny = 1;
Nz = 40;

Xrange = [-5 5]/1000; % in m
Yrange = 0;%[-0.1 0.1]/1000;
Zrange = [25 45]/1000; % in m

SimulationBox = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz);
% Hf1 = figure(1);
% set(Hf1,'name','position of detection')
% scatter3(SimulationBox.X*1e3,SimulationBox.Y*1e3,SimulationBox.Z*1e3)
% xlabel('x(mm)')
% ylabel('y(mm)')
% zlabel('z(mm)')

% calculation of the emitted Field :
[h,t] = calc_hp(Probe,SimulationBox.Points());
h = h/max(h(:));
time = t + [0:(size(h,1)-1)]/fs;
% h : each column corresponds to the field calculated for the
% SimulationBox.Points() list

% Screen Raw data Result
% Hf2 = figure(2);
% set(Hf2,'name','raw data returned By calc_hp function')
% %imagesc(1:Nx*Ny*Nz,time,log(abs(h)))
% [~,I] = sort(SimulationBox.X(:).^2 + SimulationBox.Y(:).^2 +SimulationBox.Z(:).^2);
% imagesc(1:Nx*Ny*Nz,time*1e6,log(h(:,I).^2))
% title('image in log scale')
% xlabel('Point index')
% ylabel('time (\mu s)')
% colorbar

% Screen maximum field for each position:
Field_resized = reshape(h',[Ny,Nx,Nz,length(time)]);
% first dimension now corresponds to the y values.

MaxAxis = max(Field_resized(:).^2);
[X,Y,Z] = meshgrid(SimulationBox.x*1e3,SimulationBox.y*1e3,SimulationBox.z*1e3);

xslice = [];
yslice = 0; 
zslice = [];
% 
Hf2 = figure(3);
set(Hf2,'name','field resized to simulation box')
set(Hf2,'NextPlot','replace')
% initialize frames for movie capture :


% for loop = 1:length(time)
% surfc(SimulationBox.x*1e3,SimulationBox.z*1e3,abs(squeeze(Field_resized(1,:,:,loop))).^2');
% shading interp
% %slice(X,Y,Z,abs(Field_resized(:,:,:,loop)).^2,xslice,yslice,zslice)
% %contourslice(X,Y,Z,abs(Field_resized(:,:,:,loop)).^2,cvals );
% %isosurface(X,Y,Z,abs(Field_resized(:,:,:,loop)).^2,0.1);
% title(['image in log scale at ',num2str(time(loop)),'s'])
% xlabel('x (mm)')
% ylabel('z (mm)')
% %zlabel('z (mm)')
% zlim([0 MaxAxis])
% %caxis([0 MaxAxis])
% colorbar
% drawnow
% 
% end

%screen maximum value for the field :
Field_max= reshape(max(h,[],1),[Nx,Ny,Nz]);
Hf3 = figure(3);
set(Hf3,'name','maximum field values')
%imagesc(SimulationBox.x*1e3,SimulationBox.z*1e6,squeeze(Field_max(1,:,:)));
imagesc(SimulationBox.x*1e3,SimulationBox.z*1e3,squeeze(Field_max(:,1,:))');
xlabel('x (mm)')
ylabel('z (mm)')
colorbar


xdc_free(Probe);
field_end;
