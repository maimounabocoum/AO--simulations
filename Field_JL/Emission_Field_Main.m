%%%%%%% Test simu Field 1 08_03_2013 %%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%

%%%% General Parameters %%%%%%%%%%%%%%%%%%%%%%%%
field_init(0);
savePression =1;        %Do you want to save data?
times = 0;
P0 = 40e5;              % Pression mesure au foyer en manip

%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%
R = 12.5/1000;            % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 1/1000;      % Size of math elements
f0 = 4.3e6;             % Transducer center frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
%No_elements = (R./ele_size);%.*(R./ele_size+1)       % # of elements

%%%% Field Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
Rho = 1000;             % Density
fs = 200e6;             % Sampling frequency [Hz]
attenuation = 0.6;      % in db/cm/Mhz
set_sampling(fs);
set_field('c',c);
set_field('Freq_att',attenuation*100/1e6);
set_field('att',2.6*100);
set_field('use_att',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot concave source %%%%%%%%%%%%%%%%%%%%%%%

r = R; % Plot concave source
[theta,phi] = meshgrid(linspace(0,pi,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x1 = r.*cos(theta).*cos(phi);
y1 = r.*sin(theta).*cos(phi);
z1 = r.*sin(phi);
theta = 0;
[r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x2 = r.*cos(theta).*cos(phi);
y2 = r.*sin(theta).*cos(phi);
z2 = r.*sin(phi);
theta = pi;
[r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x3 = r.*cos(theta).*cos(phi);
y3 = r.*sin(theta).*cos(phi);
z3 = r.*sin(phi);
x = [x1;x2;x3];
y = [y1;y2;y3];
z = [z1;z2;z3];

figure(1000)
subplot(2,2,1); surf(x,y,z)
axis square
view(10,100)
axis([ -R R 0 Rfocal -R R ])
axis off
subplot(2,2,2); surf(x,y,z)
axis square
view(-180,90)
axis([ -R R 0 Rfocal -R R ])
subplot(2,2,3); surf(x,y,z)
axis square
view(90,0)
axis([ -R R 0 Rfocal -R R ])
subplot(2,2,4); surf(x,y,z)
axis square
view(0,180)
axis([ -R R 0 Rfocal -R R ])
saveas(gcf,['Source_1'], 'jpg')
close(gcf)
clear x1 x2 x3 y1 y2 y3 z1 z2 z3 x y z r theta phi


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfoBeam.XDebut  = -2*R;  InfoBeam.XFin  = 2*R;  InfoBeam.pasX  = 2/1000;
InfoBeam.ZDebut  = -2*R;  InfoBeam.ZFin  = 2*R;    InfoBeam.pasZ  = 2/1000;
InfoBeam.YDebut  = 2*R;  InfoBeam.YFin  = 6*R;   InfoBeam.pasY  = 4/1000;

[points,ampl,InfoBeam] = ZoneImageRemote(InfoBeam);

[champs ,ti] = calc_hp(Th,points);
MaxChamps = max(max(champs))
champs = champs/max(max(champs));
P0 = 40e5; % Pression mesure au foyer en manip
champs = champs*P0;
champsM = max(champs,[],1);

Z = Rho*c;
C = attenuation*f0*1e-6*1e2/8.7/Rho/c^2;

Force = champsM;

Dpl = reshape(Force,InfoBeam.NbY,InfoBeam.NbX,InfoBeam.NbZ);
champF=zeros(InfoBeam.NbY,InfoBeam.NbX,InfoBeam.NbZ);
champF=Dpl;
figure;
for i=1:InfoBeam.NbY
    figure(i)
    toto = (squeeze(champF(i,:,:)))';
    imagesc(toto);  colorbar; title('champ diffracté ds l''espace');
    %colorbar
    %caxis([0 max(max(max(champF)))])
    saveas(gcf,['PressureField_1a_' num2str(i)], 'jpg')
    pause;
    close(gcf)
     %champF(:,:,i)=(champF(:,:,i)).^2; %champF(:,:,i) = champF(:,:,i)./max(max(champF(:,:,i))); 
end
for i=1:InfoBeam.NbX
    figure(i)
    toto = (squeeze(champF(:,i,:)))';
    imagesc(toto);  colorbar; title('champ diffracté ds l''espace');
    %colorbar
    %caxis([0 max(max(max(champF)))])
    saveas(gcf,['PressureField_1b_' num2str(i)], 'jpg')
    pause;
    close(gcf)
     %champF(:,:,i)=(champF(:,:,i)).^2; %champF(:,:,i) = champF(:,:,i)./max(max(champF(:,:,i))); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enregistrement des données pour la palpation simulée
name = ['PressureFiled_1.mat'];

if savePression == 1   
    fid=fopen([name '.bin'],'wb');
    fid2=fopen(['zone' name '.bin'],'wb');
    fwrite(fid,Force,'float');
    fwrite(fid2,points,'float');
    save(['Info' name],'InfoBeam');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdc_free(Th);
field_end;
fclose all;
