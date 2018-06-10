%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clearvars;
% addpath('functions')
% addpath('..\Scan routines')
% addpath('shared functions folder')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of the OP structure
                            % X : theta 
                            % Y : monotonic t variable in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% experiemental input datas :
% load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
 c = 1540 ; % sound velocity in m/s

    R_FT = MyImage.fourier(MyImage.R) ;
    FILTER = MyImage.GetFILTER(1e-3) ;
    I = MyImage.ifourier(R_FT.*FILTER) ;
    [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 1 - find x position for t = 0 :
 % Need to implement on real experiement : finding zero . Here it is not
 % necessary since t = 0 matches xsonde = 0 for all angles
 % imagesc(MyImage.theta*180/pi,xsonde, DelayLAWS) %-- view simulated delay

 X_m = (1:192)*0.2*1e-3 ;  
[theta,M0,X0,Z0] = EvalDelayLaw_shared(X_m,DelayLAWS,ActiveLIST,c); 

 Hf = figure;
Ireconstruct = Retroprojection_shared( I , X_m, z_out , theta, M0 , Hf);

% subplot(211);plot(X_m*1e3,Ireconstruct(105,:)); subplot(212);plot(-20.5 + z_out*1e3,Ireconstruct(:,96))
%  
% fwhm = FWHM(Ireconstruct(105,:),X_m*1e3)
% fwhm = FWHM(Ireconstruct(:,96),z_out*1e3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure
% imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
% colorbar
% title('simulation input phantom')
% ylim([min(z_out*1e3) max(z_out*1e3)])
% xlabel('x (mm)')
% ylabel('y (mm)')
% drawnow   



%rmpath('..\Scan routines')
%rmpath('functions')

