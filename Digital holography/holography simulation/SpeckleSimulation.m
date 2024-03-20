%% add dependencies
addpath('..\..\..\AO--commons\shared functions folder');
addpath('..\..\..\AO--commons\common subfunctions');

%% Define camera
clear all

Nx_cam = 512;   % 1024;
Ny_cam = 512;   %  620;
dpixel = 5e-6;
QE =  0.2;
Relec = 1 ; % in electrons
AD = 4.3; % electron/counts
bit = 8;
Dark = 1 ; % electron/sec/px PCO

%% US field temporal profile
N   = 2^(nextpow2( max(Nx_cam,Ny_cam) )) ;      % number of point in Fourier
F   = TF2D( N ,N, 1/(dpixel), 1/(dpixel));    % Fourier structure


%% ========= initial gaussian field in (x,y) plane + plane wave tilted

lambda0     = 780e-9;      % m
c           = 3e8;         % m/s
nu          = c/lambda0;   % Hz
h           = 6.626*1e-34; % J.s
Ephoton     = h*nu ;       % J

P0 = 6500e-6 ;
Pprobe      = (2e-4)*P0*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 );          % Power of Main pulse in W*(Scam/(1cm^2))
Pref        = P0*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 );                        % Power of Ref pulse in W*(Scam/(1cm^2))
Tint        = 100e-6;                                                               % integration time of the camera in s
ModeWidth   = 3e4 ;                                                               % reducing it the will affect final speckle size
eta         = 1 ;                                                                   % tagging efficiency (all photons are tagged when eta = 1)

% following function returns :
% E_tag   : field of tagged photon in    sqrt(W/m^2)
% E_untag : field of untagged photon in  sqrt(W/m^2)
% E_ref   : field of tagged photon in    sqrt(W/m^2)

[E0_tag,E0_untag,Eref] = initField(Pprobe,Pref,ModeWidth,eta,F);

%% intensity on the camera
I1   = Tint.*( abs(E0_untag).^2 + abs(E0_tag + 0*Eref ).^2 ); % in J/m^2
I2   = 0.5*Tint.*( abs(E0_untag).^2 + abs(E0_tag + 0*Eref ).^2 ); % in J/m^2

%% convert to J
Icam1    = I1*dpixel*dpixel ; % in J
Icam2    = I2*dpixel*dpixel ; % in J

%% ====== convert from J to Photon Number ( on each pixel) + statistical Poissoninan noise
Nphoton1        = PoissonNoise( Icam1/Ephoton ) +  PoissonNoise( 4 + 0*Icam1/Ephoton ) ; 
Nphoton2        = PoissonNoise( Icam2/Ephoton ) +  PoissonNoise( 4 + 0*Icam1/Ephoton ) ; 

%% set limit due to dark noise and saturation
Nphoton1( Nphoton1 > 2^bit )        = 2^bit ;
Nphoton1( Nphoton1 < Tint*Dark/AD ) = 0 ;
Nphoton2( Nphoton2 > 2^bit )        = 2^bit ;
Nphoton2( Nphoton2 < Tint*Dark/AD ) = 0 ;

figure(2);clf;
subplot(121)
imagesc(Nphoton1) ; colorbar ;
subplot(122)
imagesc(Nphoton2) ; colorbar ;


g = corr2(Icam1,Icam2) 
g = corr2(Nphoton1,Nphoton2) 












