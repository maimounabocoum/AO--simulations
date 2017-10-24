%%% Influence PSF on generated phantom :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 19-10-2017
clearvars ;
parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);
CurrentExperiement = CurrentExperiement.EvalPhantom();
[MyTansmission,~,~] = CurrentExperiement.ShowPhantom;
%use param.angles has an input to additionally show Radon transform
x = CurrentExperiement.MySimulationBox.x ;
z = CurrentExperiement.MySimulationBox.z ;
theta = CurrentExperiement.ScanParam ;
dx = x(2) - x(1) ;
dz = z(2) - z(1) ;

% convolution by gaussian PSF : 
% satandard deviation : f(x) = exp(-x^2/(2sigma^2))
% FWHM = 2.355 sigma
FWHMx = 1e-3 ;
FWHMz = 0.25e-3 ;
sigmax = (FWHMx/2.355) / dx  ; % in pixels
sigmaz = (FWHMz/2.355) / dz  ; % in pixels

 H = figure(101);
 set(H,'WindowStyle','docked'); 
 title('Convoluted image')
 TransmissionConv = imgaussfilt(MyTansmission,[sigmaz,sigmax]) ;
 title('idaron(radon)')
 [ImageRadon,zp] = radon(MyTansmission,theta*180/pi+ 90)  ;
 SimuImage = iradon(ImageRadon , theta*180/pi + 90) ;
 zp = (-floor(size(SimuImage,1)/2):1:floor(size(SimuImage,1)/2-1));
 x_radon = zp*dx + mean(x);
 z_radon = zp*dz + mean(z);
%x_radon*1e3,z_radon*1e3,
  imagesc(SimuImage)
%  axis([-20 20 10 35])
 %  imagesc(x*1e3,z*1e3,TransmissionConv)
 ylabel('z(mm)')
 xlabel('x(mm)')
 
%  figure(101)
%  subplot(211)
%  plot(z_radon*1e3,SimuImage(:,536))
%  subplot(212)
%  plot(x_radon*1e3,SimuImage(518,:))
% %  fwhmx = FWHM(MyTansmission(76,:),x*1e3)
% %  fwhmz = FWHM(MyTansmission(:,65),z*1e3)
%  fwhmx = FWHM(SimuImage(518,:),x_radon*1e3)
%  fwhmz = FWHM(SimuImage(:,536),z_radon*1e3)
 
 
 
 
 
 
 
 