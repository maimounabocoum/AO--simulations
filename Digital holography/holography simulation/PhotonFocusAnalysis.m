%%% analysis Victor - camera Photon MV1-D1024E-160-CL
clearvars;

Nx_cam = 1024;  % 1024;
Ny_cam = 1024;  % 620;
dpixel = 10.6e-6;
QE =  0.3 ;
Relec = 200 ; % in electrons 
AD = 48.82 ; % electron/counts
bit = 12; % 
Dark = 12500 ; % electron/sec/px PCO


%% load datas dark
Tau_exp = [11e-6,100e-6,1e-3,10e-3];
Data_dark = load('camera PhotonFocus-Victor\PhotonFocus\11us_dark_mean_std.mat');
frame0(:,:,1) = Data_dark.frame0 ;
frame1(:,:,1) = Data_dark.frame1 ;
Data_dark = load('camera PhotonFocus-Victor\PhotonFocus\100us_dark_mean_std.mat');
frame0(:,:,2) = Data_dark.frame0 ;
frame1(:,:,2) = Data_dark.frame1 ;
Data_dark = load('camera PhotonFocus-Victor\PhotonFocus\1ms_dark_mean_std.mat');
frame0(:,:,3) = Data_dark.frame0 ;
frame1(:,:,3) = Data_dark.frame1 ;
Data_dark = load('camera PhotonFocus-Victor\PhotonFocus\10ms_dark_mean_std.mat');
frame0(:,:,4) = Data_dark.frame0 ;
frame1(:,:,4) = Data_dark.frame1 ;

%% load datas Bright
Tau_exp = [11e-6,100e-6,1e-3,10e-3];
Data_bright = load('camera PhotonFocus-Victor\PhotonFocus\11us_bright_mean_std.mat');
frame0(:,:,1) = Data_bright.frame0 ;
frame1(:,:,1) = Data_bright.frame1 ;
Data_bright = load('camera PhotonFocus-Victor\PhotonFocus\100us_bright_mean_std.mat');
frame0(:,:,2) = Data_bright.frame0 ;
frame1(:,:,2) = Data_bright.frame1 ;
Data_bright = load('camera PhotonFocus-Victor\PhotonFocus\1ms_bright_mean_std.mat');
frame0(:,:,3) = Data_bright.frame0 ;
frame1(:,:,3) = Data_bright.frame1 ;
Data_bright = load('camera PhotonFocus-Victor\PhotonFocus\10ms_bright_mean_std.mat');
frame0(:,:,4) = Data_bright.frame0 ;
frame1(:,:,4) = Data_bright.frame1 ;

%% plot results
for i = 1:4
figure;
subplot(211)
imagesc(AD*frame0(:,:,i))
mu(i) = mean(mean(frame0(:,:,i)));
title(['\mu \tau_{exp}(ms) = ', num2str(Tau_exp(i)*1e3)])
cb = colorbar
ylabel(cb,'electron-counts')
subplot(212)
imagesc(AD*frame1(:,:,i))
title(['\sigma \tau_{exp}(ms) = ', num2str(Tau_exp(i)*1e3)])
std(i) = mean(mean(frame1(:,:,i)));
cb = colorbar
ylabel(cb,'electron-counts')
end

figure(2);
plot(Tau_exp,mu,'linewidth',3,'Marker','square')
