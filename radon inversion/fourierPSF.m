%%% construction of PSF through fourier FFT
clearvars

N = 2^7 ;
Fmax = 1/(0.2e-3); % 100 microns scan

FourierS = TF2D(N,Fmax);

[X,Y]   = meshgrid(FourierS.x,FourierS.y) ;
[Fx,Fy] = meshgrid(FourierS.fx,FourierS.fy) ;
THETA = atan2(Fy,Fx);

ImageDirac = zeros(size(X)) ;
% dirac
ImageDirac = exp(-X.^2/(0.2e-3)^2-Y.^2/(4e-3)^2);
% dirac
%ImageDirac(N/2+20,N/2+20) = 1;

subplot(221)
imagesc(FourierS.x*1e3,FourierS.y*1e3,ImageDirac)
title('object')
xlabel('x(mm)')
ylabel('z(mm)')

% fourier transform of pixel
ImageFourier = FourierS.fourier(ImageDirac) ;

subplot(223)
imagesc(FourierS.fx*1e-3,FourierS.fy*1e-3, real(ImageFourier))
title('fourier transform')
xlabel('fx(mm^{-1})')
ylabel('fz(mm^{-1})')

% filter in fourier domain
theta_c = 46*pi/180 ;
 FILTER = (abs(THETA-pi/2) < theta_c) | (abs(THETA+pi/2) < theta_c) ;
 FILTER(N/2+1,N/2+1) = 1 ;
% FILTER = ( sqrt(Fx.^2 + Fy.^2) < 1e3 );

subplot(224)
imagesc(FourierS.fx*1e-3,FourierS.fy*1e-3,...
        abs(ImageFourier).*FILTER)
title('filtered')
xlabel('x(mm)')
ylabel('z(mm)')

PSF = FourierS.ifourier(ImageFourier.*FILTER);

subplot(222)
imagesc(FourierS.x*1e3,FourierS.y*1e3,PSF)
title('reconstructed')
xlabel('x(mm)')
ylabel('z(mm)')

% correponding FWHM
  fwhmx = FWHM(PSF(N/2+1,:),FourierS.x*1e3)
  fwhmz = FWHM(PSF(:,N/2+1),FourierS.y*1e3)
