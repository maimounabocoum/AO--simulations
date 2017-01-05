%%%%%%% Test simu Field 1 04_03_2013 %%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Define parameters %%%%%%%%%%%%%%%%%%%

R = 12.5/1000; % Radius of the transducer [mm]
Rfocal = 40/1000; % Focal radius of the transducer [mm]
ele_size = 1/1000; % Size of math elements

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


figure(1)
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
points = [A' B' C'];
clear A B C
%C = 0:ele_size:Rfocal*2;

[hp, start_time] = calc_hp(Th,points);
start_time
D= max(max(max(hp)));

h=figure(2);
surf(hp)
shading interp
colorbar

taille_tmp = size(hp);
taille = taille_tmp(1);
h=figure(3);
% for k=1:1:taille
%     plot(hp(k,:));
%     hold on
% end
I =imagesc(hp');
%J = imrotate(I,90);
%imshow(J)
%view(180,90)
%shading interp
colorbar
PressionPascal = hp;
save('test_pressureFiled_1.mat','PressionPascal','-mat')

% h=figure(2);
% FIG = reshape(hp(1:((2*R/ele_size-1)*(2*R/ele_size-1)),1),[(2*R/ele_size-1) (2*R/ele_size-1)]);
% surf(FIG)
% %axis([0 (2*R/ele_size-1) 0 (2*R/ele_size-1) -D D])
% caxis([0 D])
% grid off
% axis off
% shading interp
% view(90,90)
% colorbar
% colorbar
% 
% h=figure(2);
% %FIG = reshape(hp(1:(R*R-2),i),[R R]);
% FIG = reshape(hp(1:((2*R/ele_size-1)*(2*R/ele_size-1)),1),[(2*R/ele_size-1) (2*R/ele_size-1)]);
% surf(FIG)
% %axis([0 (2*R/ele_size-1) 0 (2*R/ele_size-1) -D D])
% caxis([0 D])
% grid off
% axis off
% shading interp
% view(90,90)
% colorbar
% f = getframe(h);
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,(2*R/ele_size+1)) = 0;
% 
% for i=1:1:(2*R/ele_size+1)
%     h=figure(i+1);
%     FIG = reshape(hp(1:((2*R/ele_size-1)*(2*R/ele_size-1)),i),[(2*R/ele_size-1) (2*R/ele_size-1)]);
%     surf(FIG)
%     %axis([0 (2*R/ele_size-1) 0 (2*R/ele_size-1) -D D]) 
%     caxis([0 D])
%     grid off
%     axis off
%     shading interp
%     view(90,90)
%     colorbar
%     f = getframe(h);
%     im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
%     %close(h)
% end
% imwrite(im,map,'Concave_Pressure_Field.gif','DelayTime',0,'LoopCount',inf)
% close all





