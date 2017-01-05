%%%%%%% Test simu Field 1 08_03_2013 %%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%

%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%
R = 10/1000;          % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 5/1000;      % Size of math elements
f0 = 4.3e6;             % Transducer center frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
No_elements = (R./ele_size)%.*(R./ele_size+1)       % # of elements
times = 0;
P0 = 40e5;              % Pression mesure au foyer en manip

%%%% Field Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
Rho = 1000;             % Density
fs = 200e6;             % Sampling frequency [Hz]
attenuation = 0.6;      % in db/cm/Mhz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r = R; % Plot concave source
% [theta,phi] = meshgrid(linspace(0,pi,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
% x1 = r.*cos(theta).*cos(phi);
% y1 = r.*sin(theta).*cos(phi);
% z1 = r.*sin(phi);
% theta = 0;
% [r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
% x2 = r.*cos(theta).*cos(phi);
% y2 = r.*sin(theta).*cos(phi);
% z2 = r.*sin(phi);
% theta = pi;
% [r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
% x3 = r.*cos(theta).*cos(phi);
% y3 = r.*sin(theta).*cos(phi);
% z3 = r.*sin(phi);
% x = [x1;x2;x3];
% y = [y1;y2;y3];
% z = [z1;z2;z3];
% 
% figure(100)
% subplot(2,2,1); surf(x,y,z)
% axis square
% view(10,100)
% axis([ -R R 0 Rfocal -R R ])
% axis off
% subplot(2,2,2); surf(x,y,z)
% axis square
% view(-180,90)
% axis([ -R R 0 Rfocal -R R ])
% subplot(2,2,3); surf(x,y,z)
% axis square
% view(90,0)
% axis([ -R R 0 Rfocal -R R ])
% subplot(2,2,4); surf(x,y,z)
% axis square
% view(0,180)
% axis([ -R R 0 Rfocal -R R ])
% clear x1 x2 x3 y1 y2 y3 z1 z2 z3 x y z r theta phi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Emission Transducer Definition %%%%%%%%

Th = xdc_concave(R,Rfocal,ele_size) % Define concave source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Emission Signal %%%%%%%%%%%%%%%%%%%%%%%

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
RI = MakeRI_Remote(f0,fs,50);
Tpulse = length(RI)/fs;
xdc_excitation(Th,RI);
%apodisation = ones(No_elements,1);
%xdc_apodization(Th,times,apodisation');
%rhrthrth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfoBeam.XDebut  = -R;  InfoBeam.XFin  = R ;  InfoBeam.pasX  = 0.5/1000;
InfoBeam.YDebut  = -R;  InfoBeam.YFin  = R;    InfoBeam.pasY  = 0.5/1000;
InfoBeam.ZDebut  = 5*R;  InfoBeam.ZFin  = 8*R;   InfoBeam.pasZ  = 2/1000;

[points,ampl,InfoBeam] = ZoneImageRemote(InfoBeam);

size(points,1)
time = fix(clock);
datesimu = date
disp(['Biginning: ' num2str(time(4)) ' h ' num2str(time(5)) ' min ' num2str(time(6)) ' s']);
for a = 1:1:size(points,1)
    [hp ,start_time] = calc_hp(Th,points(a,:));
    champ(a,1)=max(hp);
    ti(a,1)=start_time;
end
Maxchamps = max(max(champ))
champ = champ./max(champ);
champ=champ*P0;
% champs = reshape(champ,[InfoBeam.NbX InfoBeam.NbY InfoBeam.NbZ]);

i=0;
for a = 1:1:InfoBeam.NbZ
    for b = 1:1:InfoBeam.NbX
        for c = 1:1:InfoBeam.NbY
            i=i+1;
            champs(b,c,a) = champ(i);
        end
    end
end

disp('OK2')

for i=1:1:size(champs,3)
    figure(i)
    imagesc(squeeze(champs(:,:,i)))
    %surf(champs(:,:,i))
    axis([0 InfoBeam.NbX 0 InfoBeam.NbY 0 max(max(max(champs(:,:,4:size(champs,3)))))])
    colorbar
    caxis([0 max(max(max(champs)))])
    %shading interp
    %view(0,90)
    saveas(gcf,['PressureField_00_' num2str(i)], 'jpg')
    close(gcf)
end

fgnfghfghfghfgh
Maxhp = max(max(hp))
hp = hp/max(max(hp));
hp = hp*P0;
hpM = max(hp,[],1);
Z = Rho*c;
C = attenuation*f0*1e-6*1e2/8.7/Rho/c^2;



for a = 1:1:length(A)
    for b = 1:1:length(B)
        for c = 1:1:length(C)
            points = [A(a) B(b) C(c)];
            [hp2, start_time] = calc_hp(Th,points);
%            start_time_tot(a,b,c)=start_time;
%             for k=1:1:length(hp)
%                 hp_tot(a,b,c,k)=hp(k);
%             end
            PressionPascal(a,b,c)=max(hp2);

%             figure
%             plot(hp)
        end
    end
end

%save('PressureFiled_3.mat','PressionPascal','-mat')
zrtfgzertfgzertzer
for i=1:1:length(C)
    figure(i)
    surf(PressionPascal(:,:,i))
    axis([0 length(A) 0 length(B) 0 max(max(max(PressionPascal(:,:,4:length(C)))))])
    shading interp
    view(0,90)
    saveas(gcf,['PressureField_3b_' num2str(i)], 'jpg')
    close(gcf)
end






