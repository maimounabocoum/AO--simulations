%% influence des onde US
clear all

N = 2^10 ;
Nus = 2^9 ;
wus = (15e6);
w0 = (3e8)/(790e-9);
Fmax = 1/(100e-9);
Fmaxt = 1/(10e-9);


%% US field temporal profile
P = TF1D(Nus,Fmaxt);
Profil = exp(-(P.t).^2/(10e-6)^2).*exp(-1i*wus.*P.t);
% figure; plot(P.t*1e6,real(Profil))

%% generate fourier structure
F = TF2D(N,Fmax,Fmax);
[X,Y] = meshgrid(F.x,F.y);
[FX,FY] = meshgrid(F.fx,F.fy);

%% initial gaussian field in (x,y) plane + plane wave tilted
S0 = exp(-(X.^2+Y.^2)/(2e-6)^2); 
Sref = exp(1i*4e6*X); 
%% fourier transform of initial gaussian field + phase randomization
IS0 = F.ifourier(S0);
PHASE = 1000*pi*rand(size(X)) ;
IS0 = IS0.*exp(1i*PHASE) ;
S1 = F.fourier(IS0);

Stest = F.fourier(abs(S1+Sref));




%% effect of temporal US

% for loop = 1:length(Profil)
%IS0 = IS0.*exp( 1i*real(Profil(loop))*PHASE/100 );
S1 = F.fourier(IS0);
I = abs(S1).^2;
figure(1)
subplot(121)
imagesc(F.x*1e6,F.y*1e6,I)
title(['time = ',num2str(1e6*P.t(1)),'\mu s'])
subplot(122)
hist(I)


%saveas( gcf , ['saves/image',num2str(loop),'.png'] )
% end

%% influence de la propagation
% clear all
% 
% N = 2^10 ;
% w0 = (3e8)/(1500e-9);
% Fmax = 1/(100e-9);
% F = TF2D(N,Fmax);
% [X,Y] = meshgrid(F.x,F.y);
% [FX,FY] = meshgrid(F.fx,F.fy);
% 
% % initial pulse spatial mode
% S0 = exp(-(X.^2+Y.^2)/(2e-6)^2); 
% Sref = exp(1i*4e6*X).*exp(-(X.^2+Y.^2)/(60e-6)^2); 
% IS0 = F.ifourier(S0);
% PHASE = 1*pi*rand(size(X)) ;
% IS0 = IS0.*exp(1i*PHASE) ;
% 
% IS0 = IS0.*exp(1i*PHASE);
% S1 = F.fourier(IS0);
% imagesc(F.x*1e6,F.y*1e6,abs(S1+Sref).^2)
% 
% colormap(hot)



