%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clearvars;
addpath('functions')
addpath('..\Scan routines')
addpath('shared functions folder')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of the OP structure
                            % X : theta 
                            % Y : monotonic t variable in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% experiemental input datas :
% load('experiement images - JB - test\OP0deg-2016-02-01_13-11.mat');
 c = 1540 ; % sound velocity in m/s
% MyImage = OP(data(:,:,1),X*pi/180,Y*1e-3,Param.SamplingRate*1e6,c); 

%% simulation traces 
%  load('saved images\SimulationOS.mat');
%  load('saved images\SimulationTransmissionOS.mat');
  load('saved images\SimulationOS.mat');
  load('saved images\SimulationTransmissionOS.mat');

N       = 2^12;
Lobject = 1e-3;
Fc      = 1/Lobject;  % Lobject is the size of the object to detect. Using simple model (sinc function)
                      % we set it to kc = 100/Lobject 
MyImage = MyImage.InitializeFourier(N,10*Fc);
MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
MyImage.Fmax()        % maximum frequency sampling = 1/dt
MyImage.F_R = MyImage.fourier(MyImage.R) ;
MyImage.Show_F_R();

% extract image back to initial size :
%  [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L);%MyImage.L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 1 - find x position for t = 0 :
 % Need to implement on real experiement : finding zero . Here it is not
 % necessary since t = 0 matches xsonde = 0 for all angles
 % imagesc(MyImage.theta*180/pi,xsonde, DelayLAWS) %-- view simulated delay

X_m = (0:191)*0.2*1e-3 ; 
%X_m = X_m - mean(X_m); 
% (X0,Z0) : position of inital wavefront at t=0 , M0 point on that curve
[theta,M0,X0,Z0]    = EvalDelayLawOS_shared( X_m , DelayLAWS  , ActiveLIST , c) ;

Hf = figure;
Ireconstruct = RetroprojectionOS_shared(I,X_m,ActiveLIST,z_out,theta,M0,X0,Z0,Kx,Hf);


% inversion through inverse fourier transform:




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
colorbar
title('simulation input phantom')
ylim([min(z_out*1e3) max(z_out*1e3)])
xlabel('x (mm)')
ylabel('y (mm)')
drawnow   


