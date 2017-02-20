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

%% experiemental input datas :
% load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
% c = 1540 ; % sound velocity in m/s
% MyImage = OP(data(:,:,1),X,Y,Param.SamplingRate*1e6,c); 

%% simulation traces 
load('saved images\Simulation.mat');
load('saved images\SimulationTransmission.mat');

N = 2^12;
MyImage = MyImage.InitializeFourier(N);
MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt
Lobject = 5e-3;
Fc = 100/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)
                     % we set it to kc = 100/Lobject
                   
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
% creating a fourier parameter set using the measure sample size in m

MyImage.F_R = MyImage.fourier(MyImage.R) ;

% MyImage = MyImage.PhaseCorrection(Fc);
MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening






























%% image show

% representation in polar coordinates:

%figure;
% [THETA, W] = meshgrid(MyImage.theta*pi/180,MyImage.w(N/2:end));
% [X,Y] = pol2cart(THETA, W);
% surfc(X,Y,log(abs(MyImage.F_R(N/2:end,:))))
% view(0,90)
% shading interp
% xlabel('\omega\it_{x} (\itm^{-1})')
% ylabel('\omega\it_{y} (\itm^{-1})')
% title('Fourier Transform of Radon')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 FILTER = (abs(MyImage.f).*cos(pi*MyImage.f/(2*974.0260)))'*ones(1,length(MyImage.theta));
 FILTER(abs(MyImage.f) >= Fc, :) = 0;

%MyImage.F_R(abs(MyImage.f) >= 974.0260, :)  = 0;

%  I = MyImage.ifourier(MyImage.F_R);
 I = MyImage.ifourier(MyImage.F_R.*FILTER);
  I = real(I) ;
% 
% figure;
% imagesc(real(I))

% extract image back to initial size :
 [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 xsonde = linspace(0,128*0.2e-3,128);
 xsonde = xsonde - mean(xsonde) ;
 [X,Z]= meshgrid(xsonde,z_out);
 img = zeros(size(X,1),size(X,2),'like',X);
 
  for i= 1:length(MyImage.theta)
      
        T = X.*sin( MyImage.theta(i) ) + (Z-mean(MyImage.L)).*cos( MyImage.theta(i) ) ;
        projContrib = interp1((z_out-mean(MyImage.L))',I(:,i),T(:),'linear',0);
        img = img + reshape(projContrib,length(z_out),length(xsonde)); 
       
       
       subplot(121)
%        surfc(1e3*X,1e3*Z-35,1e3*T)
%        view(0,90)
%        shading interp
       imagesc(xsonde*1e3,z_out*1e3,img)
       %imagesc(xsonde*1e3,z_out,IMG')
       title(['angle',num2str(MyImage.theta(i)*180/pi)])
       xlabel('x (mm)')
       ylabel('z (mm)')
       colorbar
       drawnow 
       
  end
  
subplot(122)
imagesc(x_phantom*1e3,z_out*1e3,MyTansmission)
title('simulation input phantom')
xlabel('x (mm)')
ylabel('y (mm)')

%rmpath('..\Scan routines')
%rmpath('functions')

