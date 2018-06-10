% plot the fourier plane scanned using OS angle/ks in polar coordinates:
clearvars;

theta = 20;
Nbx   = 0:10;

L     = 2e-3;
Lobj  = 50e-3 ;
N     = 2^10;

% frequency decimate unit
df0t = (1/Lobj) ;
ft = (-N/2:(N/2-1))*df0t;
df0s = (1/L) ;
fs = Nbx*df0s ;

H = figure('DefaultAxesFontSize',18); 
title('\le \theta \le','interpreter','latex')


Lobj = 1e-3*1000;
for i = 1:length(Nbx)
    
    [FT,THETA] = meshgrid(ft,2*pi*theta/180);
    KX = FT.*sin(THETA) + fs(i)*cos(THETA) ;
    KZ = FT.*cos(THETA) - fs(i)*sin(THETA) ;
    scatter([KX(:);-KX(:)]/1000,[KZ(:);-KZ(:)]/1000)

    axis([-4/Lobj 4/Lobj -4/Lobj 4/Lobj])

    hold on
    drawnow
end
xlabel('k_x (mm^{-1})')
ylabel('k_z (mm^{-1})')
    saveas(H,['gif folder/image',num2str(i),'.png'])

