%%%%%%% Test simu Field 1 08_03_2013 %%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%

%%%% Transducer Parameters %%%%%%%%%%%%%%%%%%%%%
R = 12.5%12.5/1000;          % Radius of the transducer [mm]
Rfocal = 40/1000;       % Focal radius of the transducer [mm]
focus = [0 0 Rfocal];   % Fixed focal point [m]
ele_size = 1/1000;      % Size of math elements
f0 = 4.3e6;             % Transducer center frequency [Hz]
c = 1540;               % Speed of sound [m/s]
lambda = c/f0;          % Wavelength
No_elements = 50;       % # of elements
times = (Rfocal./ele_size).*(Rfocal./ele_size);

%%%% Field Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
Rho = 1000;             % Density
fs = 200e6;             % Sampling frequency [Hz]
attenuation = 0.6;      % in db/cm/Mhz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = R; % Plot concave source
[theta,phi] = meshgrid(linspace(0,pi,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x1 = r.*cos(theta).*cos(phi);
y1 = r.*sin(theta).*cos(phi);
z1 = r.*sin(phi);
theta = 0;
[r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x2 = r.*cos(theta).*cos(phi);
y2 = r.*sin(theta).*cos(phi);
z2 = r.*sin(phi);
theta = pi;
[r,phi] = meshgrid(linspace(0,R,2*R/ele_size),linspace(-pi/2,pi/2,2*R/ele_size));
x3 = r.*cos(theta).*cos(phi);
y3 = r.*sin(theta).*cos(phi);
z3 = r.*sin(phi);
x = [x1;x2;x3];
y = [y1;y2;y3];
z = [z1;z2;z3];

figure(100)
subplot(2,2,1); surf(x,y,z)
axis square
view(10,100)
axis([ -R R 0 Rfocal -R R ])
axis off
subplot(2,2,2); surf(x,y,z)
axis square
view(-180,90)
axis([ -R R 0 Rfocal -R R ])
subplot(2,2,3); surf(x,y,z)
axis square
view(90,0)
axis([ -R R 0 Rfocal -R R ])
subplot(2,2,4); surf(x,y,z)
axis square
view(0,180)
axis([ -R R 0 Rfocal -R R ])
clear x1 x2 x3 y1 y2 y3 z1 z2 z3 x y z r theta phi


Th = xdc_concave(R,Rfocal,ele_size); % Define concave source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%xdc_show(Th)

A = -R:ele_size:R;
B = -R:ele_size:R;
C = 0:ele_size:2*R;


length(A)*length(B)*length(C)
length(A)
for a = 1:1:length(A)
    a
    for b = 1:1:length(B)
        
        for c = 1:1:length(C)
            points = [A(a) B(b) C(c)];
            [hp, start_time] = calc_hp(Th,points);
%            start_time_tot(a,b,c)=start_time;
%             for k=1:1:length(hp)
%                 hp_tot(a,b,c,k)=hp(k);
%             end
            PressionPascal(a,b,c)=max(hp);
%             figure
%             plot(hp)
            clear hp start_time points
        end
    end
end

%save('PressureFiled_3.mat','PressionPascal','-mat')
zrtfgzertfgzertzer
for i=1:1:length(C)
    figure(i)
    surf(PressionPascal(:,:,i))
    axis([0 length(A) 0 length(B) 0 max(max(max(PressionPascal(:,:,4:length(C)))))])
    shading interp
    view(0,90)
    saveas(gcf,['PressureField_3b_' num2str(i)], 'jpg')
    close(gcf)
end






