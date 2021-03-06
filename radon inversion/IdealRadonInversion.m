%%% inversion using ideal radon transform form original image %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
addpath('functions')
addpath('..\Scan routines')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation traces 
% load('saved images\Simulation.mat');
% load('saved images\SimulationTransmission.mat');

N = 2^12;
Lobject = 1e-3;          % object wavelength
Fc = (2*pi)/Lobject;     % Lobject is the size of the object to detect. Using simple model (sinc function)
Zref = mean(MyImage.L) ; % Z position supposed to be invariant in rotation
 
MyImage = MyImage.InitializeFourier(N,3*Fc);
%MyImage.Show_R();      % show Radon transform (ie interpolated raw data)
MyImage.Fmax()  ;       % maximum frequency sampling = 1/dt

% test dimensions:
% MyTansmission  =  0*MyTansmission;
% MyTansmission(35:55,75:95) = 1;
% before radon transform, the x and y axis must be on the same scale:
[Xphant Yphant] = meshgrid(x_phantom,z_phantom);
x_phantom = (min(x_phantom)):(z_phantom(2) - z_phantom(1)):(max(x_phantom));
[Xint Yint] = meshgrid(x_phantom,z_phantom);
MyTansmission = interp2(Xphant,Yphant,MyTansmission,Xint,Yint,'linear',0);
%%%%%%%%%%%%%% radon transform of image %%%%%%%%%%%%%%%
 [R z1]         = radon(MyTansmission,MyImage.theta*180/pi-90); %correction of default 0 definition in radon
 [iR iRfilter]  = iradon(R,MyImage.theta*180/pi-90,length(x_phantom)); 
 zR             = z1*(z_phantom(2) - z_phantom(1)) + Zref; 
     
     figure
     imagesc(MyImage.theta*180/pi,zR*1e3,R)
     xlabel('\theta(�)')
     ylabel('z(mm)')
     title('Radon Trandform of Phantom')
 
 MyRadon = interp1(zR,R,MyImage.t,'linear','extrap');
 MyRadonTF = MyImage.fourier(MyRadon) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'ram-lak';%'ram-lak' 
filt = FilterRadon(MyImage.f, MyImage.N ,FilterType , 100*Fc);

%filt = PerfectdesignFilter('hann', length(MyImage.f)/2, Fc/max(MyImage.f));
filt = filt(:);
hold on
plot(MyImage.f,filt)


 FILTER = filt*ones(1,length(MyImage.theta));
% FILTER(abs(MyImage.f) >= Fc, :) = 0;





%p = bsxfun(@times, p, H); % faster than for-loop
% I = MyImage.ifourier(MyRadonTF);
 I = MyImage.ifourier(MyRadonTF.*FILTER);

% extract image back to initial size :
 [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);
 
 figure
 imagesc(MyImage.theta*180/pi,z_out*1e3,I)
 xlabel('\theta(�)')
 ylabel('z(mm)')
 title('filter radon transform')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 xsonde = (z_out-Zref);
% xsonde = xsonde - mean(xsonde) ;
 
 % retreive t0 correction :

 % 1 - find x position for t = 0 :
 % Need to implement on real experiement : finding zero . Here it is not
 % necessary since t = 0 matches xsonde = 0 for all angles
 % imagesc(MyImage.theta*180/pi,xsonde, DelayLAWS) %-- view simulated delay

 [X,Z]= meshgrid(xsonde,z_out);
 img = zeros(size(X,1),size(X,2),'like',X);
 
  for i= 1:length(MyImage.theta)  
      % For Ideal plane Waves reconstruction
       T = X.*sin( MyImage.theta(i) ) + (Z-Zref).*cos( MyImage.theta(i) ) ;         
      % common interpolation:  
        projContrib = interp1((z_out-Zref)',I(:,i),T(:),'linear',0);
      % retroprojection:  
        img = img + reshape(projContrib,length(z_out),length(xsonde)); 
      %%% real time monitoring %%%   
       subplot(121)
       imagesc(xsonde*1e3,z_out*1e3,img)
       title(['angle',num2str(MyImage.theta(i)*180/pi)])
       xlabel('x (mm)')
       ylabel('z (mm)')
       colorbar
       drawnow        
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(121)
imagesc(xsonde*1e3,z_out*1e3,img)
title('resconstructed phantom')
xlabel('x (mm)')
xlim([min(x_phantom*1e3) max(x_phantom*1e3)])
ylim([min(z_phantom*1e3) max(z_phantom*1e3)])
ylabel('z (mm)')
colorbar


subplot(122)
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)%MyTansmission
colorbar
title('Input phantom')
ylim([min(z_out*1e3) max(z_out*1e3)])
xlabel('x (mm)')
ylabel('y (mm)')
drawnow   


