% plot the fourier plane scanned using OS angle/ks in polar coordinates:
clearvars;

theta = [-30:30];
Nbx   = [10];
L     = 10e-2;
Lobj  = 1e-3 ;
N = 2^12;

% frequency decimate unit
df0t = (1/(3*L)) ;
ft = (-N/2:(N/2-1))*df0t;
% ft = (-N/2:0)*df0t;
df0s = (1/(0.05*L)) ;
%fs = Nbx*df0s ;
fs = [1000] ;

% figure('DefaultAxesFontSize',18); 
% title('Fourier plane')

Lobj = 1e-3;
for i = 1:length(fs)
    
    [FT,THETA] = meshgrid(2*pi*ft,pi*theta/180);
    FX = FT.*sin(THETA) + fs(i)*cos(THETA) ;
    FZ = FT.*cos(THETA) - fs(i)*sin(THETA) ;
   % scatter(1e-3*FX(:),1e-3*FZ(:))
    scatter(1e-3*[FX(:);-FX(:)],1e-3*[FZ(:);-FZ(:)])%,[],repmat([0,0.45,0.74],8192,1))

    axis equal
    xlabel('f_x(mm^{-1})')
    ylabel('f_z(mm^{-1})')
    hold on
    drawnow
end

    
    xlabel('f_x(mm^{-1})')
    ylabel('f_z(mm^{-1})')
    axis([-4/Lobj 4/Lobj -4/Lobj 4/Lobj])
    
R = 4e-3*1000; %radius
S=15;   %num circ.lines
N=16;   %num ang.lines

sect_width=2*pi/N;    
offset_angle=(deg2rad(11.25)):sect_width:2*pi-sect_width+deg2rad(11.25);

%------------------
r=linspace(0,R,S+1);
w=0:.01:2*pi;

hold on
axis equal
for n=2:length(r)

      plot(real(r(n)*exp(1i*w)),imag(r(n)*exp(1i*w)),'k--')

end 

% 
for n=1:length(offset_angle)

      plot(real([0 R]*exp(1i*offset_angle(n))),imag([0 R]*exp(1i*offset_angle(n))),'k-')

end

axis(1e-3*[-(2/Lobj) (2/Lobj) -(2/Lobj) (2/Lobj)])
   set(findall(gcf,'-property','FontSize'),'FontSize',20) 


