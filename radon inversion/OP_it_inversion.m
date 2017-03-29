%%%%%%%%%%%%%%%%% using conversion interation %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
addpath('functions')
addpath('..\Scan routines')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of the OP structure
                            % X : theta 
                            % Y : monotonic t variable in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation traces 
load('saved images\Simulation.mat');
load('saved images\SimulationTransmission.mat');

%% experiemental input datas :
% load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
% c = 1540 ; % sound velocity in m/s
% MyImage = OP(data(:,:,1),X,Y,Param.SamplingRate*1e6,c); 

%% simulation traces 
load('saved images\Simulation.mat');
load('saved images\SimulationTransmission.mat');

% imported variables :
% MyTansmission : phatom to analyse

N = 2^10;
Lobject = 1e-3;
Fc = 2/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)

MyImage = MyImage.InitializeFourier(N,Fc);
%MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt

%%%%%%%%%%%%%% radon transform of image %%%%%%%%%%%%%%%
 [R z1] = radon(MyTansmission,MyImage.theta*180/pi-90); %correction of default 0 definition in radon
 zR = z1*(z_phantom(2) - z_phantom(1));   
 MyRadon = interp1(zR,R,MyImage.t,'linear',0);
 MyRadonTF = MyImage.fourier(MyRadon) ;

Original = TF2D(N,Fc);
 [Xf,Yf] = meshgrid(Original.x,Original.y) ;
 z_out = ReduceDataSize(MyImage.t,'x',MyImage.t,MyImage.L);
[X Y] = meshgrid(x_phantom,z_phantom);
 MyTansmissionInterp = interp2(X,Y, MyTansmission,Xf,Yf,'linear',0);
 MyTansmissionTF = Original.fourier(MyTansmissionInterp) ;


 MyImage.F_R = MyImage.fourier(MyImage.R);
%MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening

%% image show

% representation in polar coordinates:

figure;

subplot(2,3,1)
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
xlabel('x(mm)')
ylabel('y(mm)')

subplot(2,3,4)
imagesc(Original.fx,Original.fy,log(abs(MyTansmissionTF)))
axis(0.5*[-Fc Fc -Fc Fc])
view(0,90)
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('FT original image')

subplot(2,3,2)
imagesc(MyImage.theta*180/pi,MyImage.t*1e3+34,MyRadon)
ylim([24 44])
xlabel('\theta (°)')
ylabel('z (mm)')

subplot(2,3,5)
[THETA, F] = meshgrid(MyImage.theta-pi/2,MyImage.f(N/2:end));
[FX,FY] = pol2cart(THETA, F);
surfc(FX,FY,log(abs(MyRadonTF(N/2:end,:))))
view(0,90)
axis(0.5*[-Fc Fc -Fc Fc])
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of Th Radon')

subplot(2,3,3)
imagesc(MyImage.theta*180/pi,MyImage.t*1e3,MyImage.R)
ylim([24 44])
xlabel('\theta (°)')
ylabel('z (mm)')

subplot(2,3,6)
[THETA, F] = meshgrid(MyImage.theta-pi/2,MyImage.f(N/2:end));
[FX,FY] = pol2cart(THETA, F);
surfc(FX,FY,log(abs(MyImage.F_R(N/2:end,:))))
view(0,90)
axis(0.5*[-Fc Fc -Fc Fc])
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of measured Radon')


























%rmpath('..\Scan routines')
%rmpath('functions')

