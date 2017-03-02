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


N = 2^10;
MyImage = MyImage.InitializeFourier(N);
%MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt
Lobject = 0.8e-3;
Fc = 1/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)


Original = TF2D(N,Fc);
 [Xf,Yf] = meshgrid(Original.x,Original.y) ;
 z_out = ReduceDataSize(MyImage.t,'x',MyImage.t,MyImage.L);
[X Y] = meshgrid(x_phantom,z_out);
 MyTansmission = interp2(X,Y, MyTansmission,Xf,Yf,'linear',0);
 MyTansmissionTF = Original.fourier(MyTansmission) ;


MyImage.F_R = MyImage.fourier(MyImage.R);
%MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening

%% image show

% representation in polar coordinates:

figure;
subplot(1,2,1)
[THETA, W] = meshgrid(MyImage.theta,MyImage.w(N/2:end));
[X,Y] = pol2cart(THETA, W);
surfc(X,Y,abs(MyImage.F_R(N/2:end,:)))
axis([-Fc Fc -Fc Fc])
view(0,90)
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of Radon in polar')

subplot(1,2,2)
imagesc(Original.kx,Original.ky,abs(MyTansmissionTF))
axis([-Fc Fc -Fc Fc])

view(0,90)
shading interp
xlabel('\omega\it_{x} (\itm^{-1})')
ylabel('\omega\it_{y} (\itm^{-1})')
title('Fourier Transform of Radon in polar')
























%rmpath('..\Scan routines')
%rmpath('functions')

