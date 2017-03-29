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
Lobject = 1e-3;
Fc = 2/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)
                   % we set it to kc = 100/Lobject 
MyImage = MyImage.InitializeFourier(N,Fc);
%MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()       % maximum frequency sampling = 1/dt
                 
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
% creating a fourier parameter set using the measure sample size in m

MyImage.F_R = MyImage.fourier(MyImage.R) ;

% MyImage = MyImage.PhaseCorrection(Fc);
% MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'cosine';%'ram-lak' 

FILTER = FilterRadon(MyImage.f, MyImage.N ,FilterType , Fc);

%cos(pi*MyImage.f/(2*974.0260)))'
 FILTER = FILTER'*ones(1,length(MyImage.theta));
% FILTER(abs(MyImage.f) >= Fc, :) = 0;

%p = bsxfun(@times, p, H); % faster than for-loop
%  I = MyImage.ifourier(MyImage.F_R);
 I = MyImage.ifourier(MyImage.F_R.*FILTER);
% I = real(I) ;
% 
% figure;
% imagesc(real(I))

% extract image back to initial size :
 [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 xsonde = (1:size(DelayLAWS,1))*0.2e-3; %linspace(0,192*0.2e-3,128);
 xsonde = xsonde - mean(xsonde) ;
 
 % retreive t0 correction :
 Zref = mean(MyImage.L) ; % Z position supposed to be invariant in rotation
 % 1 - find x position for t = 0 :
 % Need to implement on real experiement : finding zero . Here it is not
 % necessary since t = 0 matches xsonde = 0 for all angles
 % imagesc(MyImage.theta*180/pi,xsonde, DelayLAWS) %-- view simulated delay

 


 [X,Z]= meshgrid(xsonde,z_out);
 img = zeros(size(X,1),size(X,2),'like',X);
 
  for i= 1:length(MyImage.theta)
      
      % For Ideal plane Waves reconstruction
       %T = X.*sin( MyImage.theta(i) ) + (Z-Zref).*cos( MyImage.theta(i) ) ;
       
      % for FIELD II reconstruction :
        T = (X-Zref*sin(MyImage.theta(i))).*sin( MyImage.theta(i) ) + (Z-Zref*cos(MyImage.theta(i))).*cos( MyImage.theta(i) ) ;
        
      % common interpolation:  
        projContrib = interp1((z_out-Zref)',I(:,i),T(:),'linear',0);
      % retroprojection:  
        img = img + reshape(projContrib,length(z_out),length(xsonde)); 
      
      %%% real time monitoring %%%   
%        subplot(121)
%        imagesc(xsonde*1e3,z_out*1e3,img)
%        title(['angle',num2str(MyImage.theta(i)*180/pi)])
%        xlabel('x (mm)')
%        ylabel('z (mm)')
%        colorbar
%        drawnow 
%        
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(121)
imagesc(xsonde*1e3,z_out*1e3,img)
title('resconstructed profile')
xlabel('x (mm)')
xlim([min(x_phantom*1e3) max(x_phantom*1e3)])
ylabel('z (mm)')
colorbar


subplot(122)
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
colorbar
title('simulation input phantom')
ylim([min(z_out*1e3) max(z_out*1e3)])
xlabel('x (mm)')
ylabel('y (mm)')
drawnow   

% subplot(223)
% title('simulation input phantom')
% 
% imvisTF = ObjectInitial.fourier(imvis);
% plot(MyImage.f,abs(imvisTF))
% xlabel('fx (m-1)')
% ylabel('y (mm)')
% XLIM = get(gca,'xlim');
% 
% subplot(224)
% ObjectInitial_FI = ObjectInitial.fourier(ObjectInitial_I);
% plot(ObjectInitial.f,sum(abs(ObjectInitial_FI),2))
% title('object fourier transform in the vertical direction')
% xlabel('fx (m-1)')
% ylabel('y (mm)')
% xlim(XLIM)


%rmpath('..\Scan routines')
%rmpath('functions')

