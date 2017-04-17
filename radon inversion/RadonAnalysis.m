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
 c = 1540 ; % sound velocity in m/s
%% simulation traces 
load('saved images\Simulation.mat');
load('saved images\SimulationTransmission.mat');

%% creation of an image using perfect radon transform :
MyPerfectImage = OP(R,MyImage.theta,zR,10e6,c) ;

N = 2^12;
Lobject = 1e-3;
Fc = 30/Lobject;    % Lobject is the size of the object to detect. Using simple model (sinc function)
MyPerfectImage = MyPerfectImage.InitializeFourier(N,Fc);
MyPerfectImage.Fmax()       % maximum frequency sampling = 1/dt
       
MyPerfectImage.Show_R();
%MyImage.Show_R();


MyPerfectImage.F_R = MyPerfectImage.fourier(MyPerfectImage.R) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'cosine';%'ram-lak' 

FILTER = FilterRadon(MyPerfectImage.f, MyPerfectImage.N ,FilterType , Fc);

%cos(pi*MyImage.f/(2*974.0260)))'
 FILTER = FILTER'*ones(1,length(MyPerfectImage.theta));
% FILTER(abs(MyImage.f) >= Fc, :) = 0;

%p = bsxfun(@times, p, H); % faster than for-loop
%   I = MyPerfectImage.ifourier(MyPerfectImage.F_R);
  I = MyPerfectImage.ifourier(MyPerfectImage.F_R.*FILTER);
% I = real(I) ;
% 
% figure;
% imagesc(real(I))

% extract image back to initial size :
 [I,z_out] = ReduceDataSize( I,'y',MyPerfectImage.t,MyPerfectImage.L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 xsonde = linspace(0,192*0.2e-3,128);
 xsonde = xsonde - mean(xsonde) ;
 Zref = mean(MyPerfectImage.L) ; % Z position supposed to be invariant in rotation

 [X,Z]= meshgrid(xsonde,z_out);
 img = zeros(size(X,1),size(X,2),'like',X);
 
  for i= 1:length(MyPerfectImage.theta)
      
      % For Ideal plane Waves reconstruction
       T = X.*sin( MyPerfectImage.theta(i) ) + (Z-Zref).*cos( MyPerfectImage.theta(i) ) ;
       
      % for FIELD II reconstruction :
      %  T = (X-Zref*sin(MyImage.theta(i))).*sin( MyImage.theta(i) ) + (Z-Zref*cos(MyImage.theta(i))).*cos( MyImage.theta(i) ) ;
        
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
% renormalization :
img = img*pi/(2*length(MyImage.theta));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% show the radon tranform of the image 
    subplot(131)
    imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
    xlabel('x(mm)')
    ylabel('z (mm)')
    xlim([min(x_phantom*1e3) max(x_phantom*1e3)])
    ylim([min(z_phantom*1e3) max(z_phantom*1e3)])

    colorbar

    subplot(132)
    imagesc(xsonde*1e3,z_out*1e3,img)
    title('reconstructed profile')
    xlabel('x (mm)')
    xlim([min(x_phantom*1e3) max(x_phantom*1e3)])
    ylim([min(z_phantom*1e3) max(z_phantom*1e3)])
    ylabel('z (mm)')
    colorbar

    subplot(133)
    [iR,H] = iradon(R,MyImage.theta*180/pi-90);
    imagesc((1:size(iR,2))*MyImage.dt*1e3,(1:size(iR,1))*MyImage.dt*1e3,iR)
    title('iverse radon transform')
    xlabel('x (mm)')
    ylabel('z (mm)')
    colorbar
