%%%%%%% Simu Field II 09_04_2013 %%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%
field_init(0);
savePression =0;        %Do you want to save data?
CurrentDirectory = cd;
%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%
R = 10/1000;            % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 1/1000;      % Size of math elements
f0 = 4.3e6;             % Transducer center frequency [Hz]
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
set_sampling(fs);
set_field('c',c);
set_field('Freq_att',attenuation*100/1e6);
set_field('att',2.6*100);
set_field('use_att',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Emission Transducer Definition %%%%%%%%

Th = xdc_concave(R,Rfocal,ele_size) % Define concave source
xdc_show(Th,'elements');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Emission Signal %%%%%%%%%%%%%%%%%%%%%%%

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
RI = MakeRI_Remote(f0,fs,50);
Tpulse = length(RI)/fs;
xdc_excitation(Th,RI);
apodisation = ones(No_elements,1);
xdc_apodization(Th,times,apodisation');
% xdc_focus(Th,times,[0 0 10]/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfoBeam.XDebut  = -0.8;  InfoBeam.XFin  = 0.8;  InfoBeam.pasX  = .2;
InfoBeam.YDebut  = -0.8;  InfoBeam.YFin  = 0.8;     InfoBeam.pasY  = .2;
InfoBeam.ZDebut  = 38;  InfoBeam.ZFin  = 42;    InfoBeam.pasZ  = 1;

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

champFtmp=zeros(InfoBeam.NbX,InfoBeam.NbY,InfoBeam.NbZ);
champFtmp=Dpl;
Plus_PetiteF = min(min(min(champFtmp)))

% Extension du FOV
InfoBeam_small = InfoBeam;
InfoBeam.XDebut = -6;
InfoBeam.XFin = 6;
InfoBeam.YDebut = -6;
InfoBeam.YFin = 6;
InfoBeam.ZDebut = 38;
InfoBeam.ZFin = 42;
Dis1X = abs(InfoBeam.XDebut - InfoBeam_small.XDebut);
Dis2X = abs(InfoBeam.XDebut - InfoBeam_small.XFin);% + InfoBeam_small.NbX;

Dis1Y = abs(InfoBeam.YDebut - InfoBeam_small.YDebut);
Dis2Y = abs(InfoBeam.YDebut - InfoBeam_small.YFin);% + InfoBeam_small.NbY;
Dis1Z = abs(InfoBeam.ZDebut - InfoBeam_small.ZDebut);
Dis2Z = abs(InfoBeam.ZDebut - InfoBeam_small.ZFin);% + InfoBeam_small.NbZ;
clear points
[points,ampl2,InfoBeam] = ZoneImageRemote2(InfoBeam);
champF=ones(InfoBeam.NbX,InfoBeam.NbY,InfoBeam.NbZ).*Plus_PetiteF ;
for xx = ((Dis1X/InfoBeam.pasX)+1):1:((Dis2X/InfoBeam.pasX+1))
    for yy = ((Dis1Y/InfoBeam.pasY)+1):1:((Dis2Y/InfoBeam.pasY+1))
        for zz = ((Dis1Z/InfoBeam.pasZ)+1):1:((Dis2Z/InfoBeam.pasZ+1))
            champF(xx,yy,zz) = champF(xx,yy,zz) + champFtmp(xx - (Dis1X/InfoBeam.pasX),yy - (Dis1Y/InfoBeam.pasY),zz - (Dis1Z/InfoBeam.pasZ));     
        end
    end
end

figure;
for i=1:10:InfoBeam.NbX
    figure(i)
    toto = (squeeze(champF(i,:,:)))';
    imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(X = ' num2str(InfoBeam.XDebut + InfoBeam.pasX*(i-1)) '))']);
    shading interp
    pause;
    close(gcf) 
end
for i=1:10:InfoBeam.NbY
    figure(i)
    toto = (squeeze(champF(:,i,:)))';
    imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(Y = ' num2str(InfoBeam.YDebut + InfoBeam.pasY*(i-1)) '))']);
    shading interp
    pause;
    close(gcf)
end
for i=1:1:InfoBeam.NbZ
    figure(i)
    toto = (squeeze(champF(:,:,i)))';
    imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(Z = ' num2str(InfoBeam.ZDebut + InfoBeam.pasZ*(i-1)) '))']);
    shading interp
    pause;
    close(gcf)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enregistrement des données pour la palpation simulée

[name,directory]=uiputfile('PressureField3D_6.mat');
if isequal(name,0) || isequal(directory,0);
else
    cd(directory);
%     
%     for i=1:InfoBeam.NbX
%         figure(i)
%         toto = (squeeze(champF(i,:,:)))';
%         imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(X = ' num2str(InfoBeam.XDebut + InfoBeam.pasX*(i-1)) '))']);
%         saveas(gcf,['PressureField_x_' num2str(i)], 'jpg')
%         close(gcf)
%     end
%     for i=1:InfoBeam.NbY
%         figure(i)
%         toto = (squeeze(champF(:,i,:)))';
%         imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(Y = ' num2str(InfoBeam.YDebut + InfoBeam.pasY*(i-1)) '))']);
%         saveas(gcf,['PressureField_y_' num2str(i)], 'jpg')
%         close(gcf) 
%     end
%     for i=1:InfoBeam.NbZ
%         figure(i)
%         toto = (squeeze(champF(:,:,i)))';
%         imagesc(toto);  colorbar; title(['champ diffracté 3D ds espace (axe(Z = ' num2str(InfoBeam.ZDebut + InfoBeam.pasZ*(i-1)) '))']);
%         shading interp
%         saveas(gcf,['PressureField_z_' num2str(i)], 'jpg')
%         close(gcf)
%     end
    
%     fid=fopen([name '.bin'],'wb');
%     fid2=fopen(['zone' name '.bin'],'wb');
%     fwrite(fid,Force,'float');
%     fwrite(fid2,points,'float');
    save([name],'champF');
    save(['zone' name],'points');
    save(['Info' name],'InfoBeam');
    save(['Info_small' name],'InfoBeam_small');
    disp('Data saved')
    cd(CurrentDirectory);
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdc_free(Th);
field_end;
fclose all;






