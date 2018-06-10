% plot the fourier plane scanned using OS angle/ks in polar coordinates:
clearvars;

theta = -20:20;
Nbx   = [5];
L     = 100e-3;
Lobj  = 100e-3 ;
N = 2^10;

% frequency decimate unit
df0t = (1/Lobj) ;
ft = (-N/2:(N/2-1))*df0t;
df0s = (1/L) ;
fs = Nbx*df0s ;

figure('DefaultAxesFontSize',18); 
title('Fourier plane')

Lobj = 50e-3;
for i = 1:length(fs)
    
    [FT,THETA] = meshgrid(ft,2*pi*theta/180);
    KX = FT.*sin(THETA) + fs(i)*cos(THETA) ;
    KZ = FT.*cos(THETA) - fs(i)*sin(THETA) ;
    scatter([KX(:);-KX(:)],[KZ(:);-KZ(:)])

    axis([-4/Lobj 4/Lobj -4/Lobj 4/Lobj])
    axis equal
    xlabel('k_x(mm^{-1})')
    ylabel('k_z(mm^{-1})')
    hold on
    drawnow
end

