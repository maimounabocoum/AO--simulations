%%
addpath('..\..\..\AO--commons\shared functions folder')
%% influence des onde US
clear all

N = 2^10 ;
Nus = 2^6 ;         % number of points in temporal profile
wus = (6e6);        % US frequency in Hz
lambda = 790e-9;
w0 = (3e8)/(lambda);% Optical frequency in Hz
Fmax = 10/(lambda);   % Echantillonage spatial
Fmaxt = 2*pi*wus/20;  % Echantillonage temporel


%% US field temporal profile
P = TF_t(Nus,Fmaxt);
Profil_US = exp(-(P.t).^2/(10e-6)^2).*exp(-1i*wus.*P.t);
figure(1); plot(P.t*1e6,real(Profil_US))

%% generate fourier structure
F = TF2D( N , N , Fmax , Fmax );
[X,Y] = meshgrid(F.x,F.z);
[FX,FY] = meshgrid(F.fx,F.fz);
SingleMode = exp(- (X.^2 + Y.^2)/(lambda)^2 );

%% fourier transform of initial gaussian field + phase randomization
IS0 = F.ifourier( sqrt(SingleMode) );
PHASE = 10000*pi*rand(size(X)) ;  % random Phase generate
Speckle = F.fourier( IS0.*exp(1i*PHASE) ) ;

figure(2)
imagesc(F.x/lambda,F.z/lambda,SingleMode);
xlabel('x(\mu m)')
ylabel('z(\mu m)')
colormap(hot)

figure(3)
imagesc(F.x/lambda,F.z/lambda,abs(Speckle).^2);
xlabel('x(\mu m)')
ylabel('z(\mu m)')
colormap(hot)


%% effect of temporal US

eta = 0.1; % tagging efficiency
Tagg_efficiency = sqrt(eta)*abs(Profil_US);

 for loop = 1:length(Profil_US)
     
   
IS0 = IS0.*exp( 1i*real(Profil_US(loop))*PHASE );

Speckle_tagged = Tagg_efficiency(loop)*F.fourier(IS0);

figure(4)
imagesc( F.x/lambda, F.z/lambda, abs(Speckle + Speckle_tagged).^2)
title(['time = ',num2str(1e6*P.t(1)),'\mu s'])
colormap(hot)
drawnow 

saveas( gcf , ['Q:\datas\simulated datas\images vrac\image',num2str(loop),'.png'] )
end
 
%
% end

%% influence de la propagation
% clear all
% 
% N = 2^10 ;
% w0 = (3e8)/(1500e-9);
% Fmax = 1/(100e-9);
% F = TF2D(N,Fmax);
% [X,Y] = meshgrid(F.x,F.z);
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
% imagesc(F.x*1e6,F.z*1e6,abs(S1+Sref).^2)
% 




