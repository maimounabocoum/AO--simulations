%%=============  add path 
Include;

%% ============= simulation parameters ===%%

Nx_cam = 2^10;   % 1024;
Ny_cam = 2^10;   % 620;
dpixel = 5.3e-6;
QE =  0.3 ;
Relec = 200 ;   % in electrons 
AD = 48.82 ;    % electron/counts
bit = 12;       % 
Dark = 12500 ;  % electron/sec/px PCO

% frequency = 0 corresponds to point of coordinate N/2+1
F       = TF2D( Nx_cam , Ny_cam , 1/(dpixel), 1/(dpixel));

param1 = (10e-6)*abs( cos(1:100 + rand(1,100)) ) ; % ref W
param2 = (1e-6)*abs( rand(1,100) ) ; % main W


%% ==== simulate noise ==== %%
    s           = zeros(1,length(param1));
    s_filtered  = zeros(1,length(param1));
    s_fft       = zeros(1,length(param1));
    
for i_loop = 1:length(param1)

lambda0     = 780e-9;      % m
c           = 3e8;         % m/s
nu          = c/lambda0;   % Hz
h           = 6.626*1e-34; % J.s
Ephoton     = h*nu ;       % J
P0          = param2(i_loop)*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 );           % Power of Main pulse in W*(Scam/(1cm^2))
Pref        = param1(i_loop)*( (Nx_cam*Ny_cam*dpixel^2)/ (1e-2)^2 );  % Power of Ref pulse in W*(Scam/(1cm^2))
Tint        = 100e-6;       % integration time of the camera in s
ModeWidth   = 50 ;          % reducing it the will affect final speckle size
eta         = 1 ;           % tagging efficiency (all photons are tagged when eta = 1)


[E0_tag,E0_untag,Eref] = initField(P0,Pref,ModeWidth,eta,F);

Icam = Tint.*( abs(E0_untag).^2 + abs(E0_tag + Eref ).^2 )*(F.dx*F.dz); % in J
Nphoton        = PoissonNoise( Icam/Ephoton );            % in photons
Nelectron    =  QE*Nphoton ;                              % in electrons
Ncount     = Nelectron/AD ;                               % in counts

Ncount_fft = F.fourier( Ncount ) ;

        % selection of ROI FFT:
        nx0  = 498 ; % 440-440 
        ny0  = 628;  % 650-440
        Nx   = 30 ; % 140-140
        Ny   = 30; % 165-140
        
        % GetFilter ( nx0 , ny0, dimension filter x , dimension filter y, image for size )
        Filt = GetFilter(nx0,ny0,Nx,Ny, Ncount_fft );

Ncount_filtered = F.ifourier( Ncount_fft.*Filt ) ;   

        s(i_loop) = sum(Ncount(:)) ;
        s_fft(i_loop) = sum( abs( Ncount_fft(:).*Filt(:) ).^2 ) ;
        s_filtered(i_loop) = sum( abs(Ncount_filtered(:)).^2 ) ;
        
        
figure(2);
subplot(121)
imagesc(1e3*F.x,1e3*F.z, Ncount )
title('tagged')
cb = colorbar; ylabel(cb,'counts');

subplot(122)
imagesc(log( abs( Ncount_fft ))) % 
title('ref')
cb = colorbar ; ylabel(cb,'\mu W') ;
        rectangle('Position',[ nx0 , ny0 , Nx , Ny ],...
         'LineWidth',1,'LineStyle','--')



end




figure
% plot(param1/max(param1)); hold on 
% plot(s/max(s)); hold on 
plot(s_fft./param1/max(s_fft./param1)); hold on 
plot(param2/max(param2)); hold on 







