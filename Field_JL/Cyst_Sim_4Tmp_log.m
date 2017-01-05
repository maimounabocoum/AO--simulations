%%%%%%% Simu Field II 04_04_2013 %%%%%%%%%%%%%%%%

clear all
%close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%
field_init(0);
savePression =0;        %Do you want to save data?
CurrentDirectory = cd;
%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%
R = 3.5/1000;            % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 0.5/1000;      % Size of math elements
f0 = 5.5e6;             % Transducer center frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
No_elements = 1;        % # of elements
times = 0;
P0 = 40e5;              % Pression mesure au foyer en manip

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

N_elements = 1;%20;% 1;
height = 0.25/1000;
width = 2.5/1000;
kerf = height/4;
Rfocus = 5/1000;
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
no_sub_x = 1;
no_sub_y = 10;

%Th = xdc_focused_array(N_elements, height, width, kerf, Rfocus, no_sub_x, no_sub_y, focus);
Th = xdc_concave(R,Rfocal,ele_size) % Define concave source
xdc_show(Th,'elements');
data = xdc_get(Th,'rect');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Emission Signal %%%%%%%%%%%%%%%%%%%%%%%

impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse = impulse.*hanning(max(size(impulse)))';
xdc_impulse (Th, impulse);
RI = MakeRI_Remote(f0,fs,50);
Tpulse = length(RI)/fs;
xdc_excitation(Th,RI);
apodisation = ones(N_elements,1);
xdc_apodization(Th,times,apodisation');
% xdc_focus(Th,times,[0 0 10]/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfoBeam.XDebut  = 0;  InfoBeam.XFin  = 0;  InfoBeam.pasX  = .2;
InfoBeam.YDebut  = -10;  InfoBeam.YFin  = 10;     InfoBeam.pasY  = .25;
% InfoBeam.XDebut  = -10;  InfoBeam.XFin  = 10;  InfoBeam.pasX  = .5;
% InfoBeam.YDebut  = 0;  InfoBeam.YFin  = 0;     InfoBeam.pasY  = .5;
InfoBeam.ZDebut  = 30;  InfoBeam.ZFin  = 50;    InfoBeam.pasZ  = .25;

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
champF=(Dpl);
%champF=20*log10(Dpl./Dpl(1,1+10/0.5,80));
%champF=20*log10(Dpl./max(max(Dpl)));
%figure;

for i=1:InfoBeam.NbX
    figure;
    toto = (squeeze(champF(i,:,:)))';
    I = imagesc(toto);  colorbar; title(['champ diffracté 2D log ds espace (axe(X = ' num2str(InfoBeam.XDebut + InfoBeam.pasX*(i-1)) '))']);
    shading interp
    %pause;
    %close(gcf) 
end
% for i=1:InfoBeam.NbY
%     figure;
%     toto = (squeeze(champF(:,i,:)))';
%     imagesc(toto);  colorbar; title(['champ diffracté ds espace (axe(Y = ' num2str(InfoBeam.YDebut + InfoBeam.pasY*(i-1)) '))']);
%     shading interp
% %     pause;
% %     close(gcf)
% end
% for i=1:InfoBeam.NbZ
%     figure(i)
%     toto = (squeeze(champF(:,:,i)))';
%     imagesc(toto);  colorbar; title(['champ diffracté ds espace (axe(Z = ' num2str(InfoBeam.ZDebut + InfoBeam.pasZ*(i-1)) '))']);
%     shading interp
%     pause;
%     close(gcf)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enregistrement des données pour la palpation simulée

[name,directory]=uiputfile('PressureField_5.mat');
if isequal(name,0) || isequal(directory,0);
else
    cd(directory);
    
    for i=1:InfoBeam.NbX
        figure(i)
        toto = (squeeze(champF(i,:,:)))';
        imagesc(toto);  colorbar; title(['champ diffracté 2D log ds espace (axe(X = ' num2str(InfoBeam.XDebut + InfoBeam.pasX*(i-1)) '))']);
        saveas(gcf,['PressureField2D_log_x_' num2str(i)], 'jpg')
        close(gcf)
    end
%     for i=1:InfoBeam.NbY
%         figure(i)
%         toto = (squeeze(champF(:,i,:)))';
%         imagesc(toto);  colorbar; title(['champ diffracté ds espace (axe(Y = ' num2str(InfoBeam.YDebut + InfoBeam.pasY*(i-1)) '))']);
%         saveas(gcf,['PressureField_y_' num2str(i)], 'jpg')
%         close(gcf) 
%     end
%     for i=1:InfoBeam.NbZ
%         figure(i)
%         toto = (squeeze(champF(:,:,i)))';
%         imagesc(toto);  colorbar; title(['champ diffracté ds espace (axe(Z = ' num2str(InfoBeam.ZDebut + InfoBeam.pasZ*(i-1)) '))']);
%         shading interp
%         saveas(gcf,['PressureField_z_' num2str(i)], 'jpg')
%         close(gcf)
%     end
     
%     fid=fopen([name '.bin'],'wb');
%     fid2=fopen(['zone' name '.bin'],'wb');
%     fwrite(fid,Force,'float');
%     fwrite(fid2,points,'float');
%     save(['Info' name],'InfoBeam');
%     disp('Data saved')
%     cd(CurrentDirectory);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdc_free(Th);
field_end;
fclose all;
