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

%MyImage = MyImage.PhaseCorrection(Fc);
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

% filtered inverse fourier transform :
% defninition of the filter :
% LP=cos(pi*f/(2*BP));
% LP(abs(f)>BP)=0;

% FILTER = (abs(MyImage.f).*cos(pi*MyImage.f/(2*974.0260)))'*ones(1,length(MyImage.theta));
% FILTER(abs(MyImage.f) >= Fc, :) = 0;

%MyImage.F_R(abs(MyImage.f) >= 974.0260, :)  = 0;

I = MyImage.ifourier(MyImage.F_R);
% I = MyImage.ifourier(MyImage.F_R.*FILTER);

%% reconstruction BOX initialization :
 xsonde = linspace(0,192*0.2e-3,193);
 xsonde = xsonde - mean(xsonde) ;
 [x,z]= meshgrid(xsonde,MyImage.z);
 img = zeros(size(x,1),size(x,2),'like',x);
 
  for i=1:length(MyImage.theta)
      
        t = x.*sin( MyImage.theta(i) ) + z.*cos( MyImage.theta(i) ) ;
        projContrib = interp1(MyImage.t',I(:,i),t(:),'linear',0);
        img = img + reshape(projContrib,length(MyImage.z),length(xsonde)); 
       
       %[IMG,z_out] = ReduceDataSize(img ,'y',MyImage.t,MyImage.L);
       subplot(121)
       imagesc(xsonde*1e3,MyImage.z*1e3,img)
       %imagesc(xsonde*1e3,z_out,IMG')
       title(['angle',num2str(MyImage.theta(i)*180/pi)])
       xlabel('x (mm)')
       ylabel('z (mm)')
       colorbar
       drawnow 
       
  end
  
subplot(122)
imagesc(x_phantom*1e3,y_phantom*1e3,MyTansmission)
title('simulation input phantom')
xlabel('x (mm)')
ylabel('y (mm)')

%rmpath('..\Scan routines')
%rmpath('functions')

