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
Rho = 1000;
attenuation = 0.6;      % in db/cm/Mhz

N_elements = 1;
focus = [0 0 40]/1000;
R_elem = 2/1000;
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

[phantom_positions, phantom_amplitudes] = phantom(2000);

%%%% ROT et trans
% Rotate = [cos(atan(-0.1388)) 0 sin(atan(-0.1388)) ; 0 1 0 ; -sin(atan(-0.1388)) 0 cos(atan(-0.1388))];
% Translate = [-5.5/1000 0 0.3831/1000];
% for i=1:1:size(phantom_positions,1)
%     phantom_positions2(i,:) = phantom_positions(i,:)*Rotate + Translate;
% end

% [v,t]= calc_scat_multi(Th,Th,[0 0 20]/1000,1);
%[v,t]= calc_scat_multi(Th,Th,phantom_positions,phantom_amplitudes);
[v,t]= calc_scat_multi(Th,Th,phantom_positions,phantom_amplitudes);

% min_sample = t*fs;
% 
% rf_env=abs(hilbert([zeros(round(t*fs-min_sample),1); v(:)]));
% env(1:size(rf_env,1)) = rf_env;


%subplot(211)
figure
[N,M]=size(v);
v=v/max(max(v));
for i=1:N_elements
    plot(2*(0:N-1)/fs+t,v(:,i)+i), hold on
end
hold off
title('Individual traces')
xlabel('Time [s]')
ylabel('Normalized response')
axis([t t+2*N/fs 0 M+1])
%subplot(212)
figure
if size(v,2)>1
    plot(2*(0:N-1)/fs+t,sum(v'))
else
    plot(2*(0:N-1)/fs+t,v(:,i)+i)
end
title('Summed response')
xlabel('Time [s]')
ylabel('Normalized response')
axis([t t+2*N/fs 0 M+1])





















impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
%xdc_impulse (Th, impulse);
figure; plot(impulse)

excitation = sin(2*pi*f0*(0:1/fs:2/f0));
%xdc_excitation (Th, impulse);
figure; plot(excitation)

%[phantom_positions, phantom_amplitudes] = phantom(2000);

%%%% ROT et trans


 %[v,t]= calc_scat_multi(Th,Th,[0 0 20]/1000,1);
%[v,t]= calc_scat_multi(Th,Th,phantom_positions,phantom_amplitudes);

%min_sample = t*fs;

%rf_env=abs(hilbert([zeros(round(t*fs-min_sample),1); v(:)]));
%env(1:size(rf_env,1)) = rf_env;


% %subplot(211)
% figure
% [N,M]=size(v);
% v=v/max(max(v));
% for i=1:N_elements
%     plot((0:N-1)/fs+t,v(:,i)+i), hold on
% end
% hold off
% title('Individual traces')
% xlabel('Time [s]')
% ylabel('Normalized response')
% %subplot(212)
% figure
% if size(v,2)>1
%     plot((0:N-1)/fs+t,sum(v'))
% else
%     plot((0:N-1)/fs+t,v(:,i)+i)
% end
% title('Summed response')
% xlabel('Time [s]')
% ylabel('Normalized response')

times = 0;
%N_elements = 1;
impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
RI = MakeRI_Remote(f0,fs,50);
Tpulse = length(RI)/fs;
xdc_excitation(Th,RI);
%apodisation = ones(N_elements,1);
%xdc_apodization(Th,times,apodisation');

InfoBeam.XDebut  = 0;  InfoBeam.XFin  = 0;  InfoBeam.pasX  = .2;
InfoBeam.YDebut  = -2;  InfoBeam.YFin  = 2;     InfoBeam.pasY  = .05;
InfoBeam.ZDebut  = 30;  InfoBeam.ZFin  = 50;    InfoBeam.pasZ  = 0.5;

[points,ampl,InfoBeam] = ZoneImageRemote2(InfoBeam);

[champs ,ti] = calc_hp(Th,points);
P0 = 40e5;
y = abs(fft(champs));
champs = champs*P0/max(max(y));
champsM = max(y,[],1);

Z = Rho*c;
C = attenuation*f0*1e-6*1e2/8.7/Rho/c^2;

Force = champsM;



counter = 0;
for xx = 1:1:InfoBeam.NbX;
    for yy = 1:1:InfoBeam.NbY;
        for zz = 1:1:InfoBeam.NbZ;
            counter = counter + 1;
            Dpl(xx,yy,zz)=Force(counter);
        end
    end
end

champF=zeros(InfoBeam.NbX,InfoBeam.NbY,InfoBeam.NbZ);
%champF=Dpl;
Dpl=20*log10(Dpl./max(max(max(Dpl))));



champF=20*log10(Dpl./max(max(Dpl)));

figure;

for i=1:InfoBeam.NbX
    figure(i)
    toto = (squeeze(champF(i,:,:)))';
    imagesc(toto);  colorbar; title(['champ diffracté 2D log ds espace (axe(X = ' num2str(InfoBeam.XDebut + InfoBeam.pasX*(i-1)) '))']);
    shading interp
    %pause;
    %close(gcf) 
end




% % for i=round(size(Dpl,1)/2):1:round(size(Dpl,1)/2)
% %     figure;imagesc(squeeze(Dpl(i,:,:)))
% % end
% for i=round(size(Dpl,2)/2):1:round(size(Dpl,2)/2)
%     figure;imagesc(squeeze(Dpl(:,i,:)))
% end
% % for i=round(size(Dpl,3)/2):1:round(size(Dpl,3)/2)
% %     figure;imagesc(squeeze(Dpl(:,:,i)))
% % end 
%     figure;plot(squeeze(Dpl(round(size(Dpl,1)/2),round(size(Dpl,2)/2),:)))

