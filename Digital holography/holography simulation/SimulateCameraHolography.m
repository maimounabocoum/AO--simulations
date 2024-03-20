addpath('..\..\..\AO--commons\shared functions folder');
addpath('..\..\..\AO--commons\common subfunctions');

%% addpath
%Include;
clearvars
% clearvars  -except SIG
isIm   = 1;     % screen out images jtest

% here : param is the Ref intensity W/cm^2
Nav = 50 ;
 param = repmat(linspace(1e-6,60e-4,20),Nav,1) ;
 param = param(:);
 
% 1024px: SIG  9.1778e-09 , NOISE 5.9459e-11  = > S/N : 154.3550
% 2048px: SIG 1.4584e-07 , NOISE 9.3700e-10 => S/N : 155.6457

% nameCamera = 'PCO.edge';
nameCamera = 'Ximea';

% turn on(=1)/off(=0) noises
IsPoissonNoise      = 1 ;
ReadoutNoise        = 0 ;
DarkNoise           = 0 ;
DigitalizedNoise    = 0 ;
MinusBG             = 1 ; % remove backgroung

for i_loop = 1%:length(param)


%%=========================== define camera screen ==============
 switch nameCamera
     case 'PCO.edge'
Nx_cam = 512;   % 1024;
Ny_cam = 512;   %  620;
dpixel = 5e-6;
QE =  0.2;
Relec = 1 ; % in electrons
AD = 4.3; % electron/counts
bit = 16;
Dark = 1 ; % electron/sec/px PCO
     case 'Ximea'
Nx_cam = 512;  % 1024;
Ny_cam = 512;   %  620;
dpixel = 5e-6;
QE =  0.2;
Relec = 1 ; % in electrons
AD = 4.3; % electron/counts
bit = 16;
Dark = 1 ; % electron/sec/px PCO

    case 'PhotonFocus'
Nx_cam = 512;  % 1024;
Ny_cam = 512;  % 620;
dpixel = 10.6e-6;
QE =  0.3 ;
Relec = 200 ; % in electrons 
AD = 48.82 ; % electron/counts
bit = 12; % 
Dark = 12500 ; % electron/sec/px PCO

    case 'manual'
Nx_cam = 2048;  % 1024;
Ny_cam = 2048;  % 620;
dpixel = 5e-6;
QE =  0.2 ;
Relec = 5 ; % in electrons      
AD = 4 ; % electron/counts
bit = 16;
Dark = 10 ; % electron/sec/px PCO



end


np = 1 ;                                            % np = number of points / pixels
N = 2^(nextpow2( np*max(Nx_cam,Ny_cam) ));          % number of point in Fourier
F       = TF2D( N ,N, np/(dpixel), np/(dpixel));    % Fourier structure


% zero of camera : inf( length(.)/2 + 1 )
% renetering pixels to zero:
% Nx_center = floor(Nx_cam/2) + 1 ;
% Ny_center = floor(Ny_cam/2) + 1 ;

x_cam =  F.x( ( -(floor(Nx_cam/2)):(floor(Nx_cam/2)-1) )*np + (N/2+1) ) ;
y_cam =  F.z( ( -(floor(Ny_cam/2)):(floor(Ny_cam/2)-1) )*np + (N/2+1) ) ;
% 
% [Xc,Yc] = meshgrid(x_cam,y_cam);

%% =========================== generate fourier structure


clearvars N Fmax

%% ========= initial gaussian field in (x,y) plane + plane wave tilted

lambda0     = 800e-9;      % m
c           = 3e8;         % m/s
nu          = c/lambda0;   % Hz
h           = 6.626*1e-34; % J.s
Ephoton     = h*nu ;       % J

P0          = (2e-4)*param(i_loop)*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 );        % Power of Main pulse in W*(Scam/(1cm^2))
Pref        = param(i_loop)*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 ); % Power of Ref pulse in W*(Scam/(1cm^2))
Tint        = 100e-6;       % integration time of the camera in s
ModeWidth   = 0.2e4 ;         % reducing it the will affect final speckle size
eta         = 0.01 ;        % tagging efficiency (all photons are tagged when eta = 1)

% following function returns :
% E_tag   : field of tagged photon in sqrt(W/m^2)
% E_untag : field of untagged photon in sqrt(W/m^2)
% E_ref   : field of tagged photon in sqrt(W/m^2)
% we avec the relation doublesum( abs(E_tag).^2 + abs(E_untag).^2 ) = P0
% important remark : energy is normalized over the fourier domain !!
% test : trapz( F.x, trapz( F.y, abs(E0_tag).^2 ) )

[E0_tag,E0_untag,Eref] = initField(P0,Pref,ModeWidth,eta,F);

%   figure; imagesc(1e3*F.x,1e3*F.z,abs(E0_tag).^2)
%   xlabel('x(mm)')
%   ylabel('x(mm)')
%   cb = colorbar
%   xlabel(cb,'W/m^2')
%   title(['Tagged light on camera. P_0 = ',num2str(P0*1e9),'nW'])

%% ===== intensity on the camera
I1   = Tint.*( abs(E0_untag).^2 + abs(E0_tag + Eref ).^2 ); % in J/m^2
I_bg =  Tint.*( abs(Eref).^2 ); % in J/m^2

% std = sqrt(var) = sigma
%  I1 = PoissonNoise((F.dx*F.dy)*I/Ephoton); % in photon
%  I1 = I1*Ephoton/(F.dx*F.dy) ;  % convert back to photon count to J/m^2
%   mean(I1(:))
%  sqrt(var(I1(:)))
% so far, I1 is interpolated on F.x , F.y , therefore its size is N^2

%% ====================== Visualization of CCD camera
% Icam : Energy on each pixels in J
Icam    = SumCamera( x_cam , y_cam  , F.x , F.z , I1 , np ); % in J
Icam_bg = SumCamera( x_cam , y_cam  , F.x , F.z , I_bg , np ); % in J 

% size of Icam is still Nx_cam * Ny_cam
Icam_tagged = SumCamera( x_cam , y_cam  , F.x , F.z , Tint*abs(E0_tag).*abs(E0_tag) , np ); % in J


% add Poisson shot noise to the image
if IsPoissonNoise == 0
 Nphoton        =  Icam/Ephoton ;
 Nphoton_bg     =  Icam_bg/Ephoton ;
else
 Nphoton        = PoissonNoise( Icam/Ephoton );
 Nphoton_bg     = PoissonNoise( Icam_bg/Ephoton );
end

%  mean(Nphoton(:))
%  sqrt(var(1*Nphoton(:)))
 

disp(['Mean Obj photons/pixels = ',num2str(mean(Icam_tagged(:)./Ephoton(:)))]);

%% effet of camera 
% Quantum efficiency model by binomial law
 Nelectron    =  QE*Nphoton ; 
 Nelectron_bg =  QE*Nphoton_bg ; 
% Nelectron = PoissonNoise(QE.*Nphoton);


% readout noise
if ReadoutNoise == 1
 Nelectron      = Nelectron + randn(size(Nelectron)).*Relec;
 Nelectron_bg   = Nelectron_bg + randn(size(Nelectron_bg)).*Relec;
end

 % Dark Current% electron/sec/px DMK
 if DarkNoise == 1
 Nelectron      = Nelectron + floor( Dark*Tint + randn(size(Nelectron)).*sqrt(Dark*Tint) );
 Nelectron_bg   = Nelectron_bg + floor( Dark*Tint + randn(size(Nelectron_bg)).*sqrt(Dark*Tint) );

 end

% Nelectron(  Nelectron(:) < 0 ) = 0; 

 
% mean(Nelectron(:)) 
% sqrt(var(Nelectron(:)))

% A/D conversion factor in e/count unit
 if DigitalizedNoise == 0
  Ncount = Nelectron/AD ;
  Ncount_bg = Nelectron_bg/AD ;
 else
  Ncount = floor(Nelectron/AD) ;
  Ncount_bg = floor(Nelectron_bg/AD) ;
 end
 

% digitalization step
Ncount = min(Ncount,2^bit);
Ncount_bg = min(Ncount_bg,2^bit);


%% Fourier transform analysis
Ng = 2^( nextpow2( max(Nx_cam,Ny_cam) ) );
G = TF2D( Ng , Ng ,1/dpixel,1/dpixel);

% usefull if image pixel dimension is not of size 2^n
%  Ncount = PadImage(Ncount,G.Nx,G.Nz);
% 
%  Icam_tagged    = PadImage(Icam_tagged,G.Nx,G.Nz);
 Ntagged        = Icam_tagged/Ephoton;
 
% sumNcount = trapz(G.x,trapz(G.y,abs(Ncount).^2))
 
% Icam_tagged = interp2(Xc,Yc,Icam_tagged,Xg,Yg,'nearest',0);
Ncom_fft        = G.fourier(Ncount);
Ncom_fft_bg     = G.fourier(Ncount_bg);

% sumNfft = trapz(G.fx,trapz(G.fy,abs(Ncom_fft).^2))
%Ntagged_fft = G.fourier(Icam_tagged/Ephoton);

% mean( abs(Ncom_fft(:)) )

fx_c0 = 1.5*15920;
fz_c0 = 0;
fx_c1 = 1.5*15920;
fz_c1 = 1.5*15920;
fr = 4500 ; % 2400

[FX,FZ]     = meshgrid(G.fx,G.fz);

Filter0     = ((FX-fx_c0).^2 + (FZ-fz_c0).^2 <= (fr)^2);
Filter1     = ((FX-fx_c1).^2 + (FZ-fz_c1).^2 <= (fr)^2);

Ncom        = G.ifourier( Ncom_fft.*Filter0 );
Ncom_bg     = G.ifourier( Ncom_fft_bg.*Filter0 );
%Ntagged    = G.ifourier(Ntagged_fft.*Filter1);
% Ncom = 2*Ncom ; % normalize positive portion of fourier domaine





if isIm == 1
figure(1);
subplot(222) ; imagesc(G.x*1e3,G.z*1e3,Ncount); cb = colorbar ; ylabel(cb,'Photo-electron count')
title('Nphoton on camera')

% mean(Ncount(:))
% sqrt(var(Ncount(:)))

figure(1);
subplot(224) ; imagesc(G.fx,G.fz,log( abs(Ncom_fft))); cb = colorbar ; ylabel(cb,'TF-log')
%subplot(224) ; imagesc(G.fx,G.fy,UnwrapPHASE(angle(Ncom_fft),G.N/2,510)); cb = colorbar ; ylabel(cb,'TF-log')
title('FFT Photo-electron')
axis([-50000 50000 -50000 50000])

theta = 0 : (2 * pi / 10000) : (2 * pi);
pline_x0 = fr * cos(theta) + fx_c0;
pline_z0 = fr * sin(theta) + fz_c0;
pline_x1 = fr * cos(theta) + fx_c1;
pline_z1 = fr * sin(theta) + fz_c1;
hold on;
plot(pline_x0, pline_z0, 'r-', 'LineWidth', 3);
plot(pline_x1, pline_z1, 'r-', 'LineWidth', 3);

figure(1);
subplot(221) ;imagesc(G.x*1e3,G.z*1e3,abs(Ntagged)) ; colorbar
axis(1e3*[min(x_cam) max(x_cam) min(y_cam) max(y_cam)])
title('Nphoton tagged')

 
figure(1);
subplot(223) ; imagesc(G.x*1e3,G.z*1e3,(abs(Ncom)-abs(Ncom_bg))*AD/QE);  cb = colorbar ;
axis(1e3*[min(x_cam) max(x_cam) min(y_cam) max(y_cam)])
title('inverse FFT - Nphoton Tagged')

end

%% signal and noise
SIG(i_loop) =  trapz(G.fx,trapz(G.fz,abs(Ncom_fft.*Filter0).^2));

NOISE(i_loop) = trapz(G.fx,trapz(G.fz,abs(Ncom_fft.*Filter1).^2)); 

end

 
for i = 1:20
   mu(i) = mean( SIG( ((1+Nav*(i-1)):(Nav*(i)) ) ) );
   sigma(i) = std( SIG( ((1+Nav*(i-1)):(Nav*(i)) ) ) ,1);
   Param(i) = mean( param( ((1+Nav*(i-1)):(Nav*(i)) ) ) );
end
 

figure(3); clf
plot(SIG)

figure(4); clf
plot(mu,mu./sigma)

