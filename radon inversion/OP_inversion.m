%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% load experiemental data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clearvars;
addpath('..\..\AO--commons\shared functions folder')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of the OP structure
                            % X : theta 
                            % Y : monotonic t variable in points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% experiemental input datas :
%% ================ load datas ==============
 clearvars

%%%%%%%%%%%% target folder to save simulated data %%%%%%%%%%
SimuPathFolder = 'Q:\datas\simulated datas';
SubFolderName = generateSubFolderName(SimuPathFolder) ;

%%
%%%%%%%%% loading OP image and experiment parameters %%%%%

 load([SubFolderName,'\1mmInclusions_type_OP_17h30_54.mat'])
 c = param.c ; % sound velocity in m/s

 % low-pass filtering of acquired data
   % dimension of data in OP structure is 1024 along t , number of angles
   % along x
    R_FT = MyImage.fourier(MyImage.R) ;     % fourier transform along radius
    FILTER = MyImage.GetFILTER(1e-3) ;      % enter cut-off length
    I = MyImage.ifourier(R_FT.*FILTER) ;    % ifourier transform along radius
    [I,z_out] = ReduceDataSize( I,'y',MyImage.t,MyImage.L); % size to original data size to actual windows of interest

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% reconstruction BOX initialization (retroprojection):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % 1 - find x position for t = 0 :
 % Need to implement on real experiement : finding zero . Here it is not
 % necessary since t = 0 matches xsonde = 0 for all angles
 % imagesc(MyImage.theta*180/pi,xsonde, DelayLAWS) %-- view simulated delay

 %%
 
 % building the X axis of the reconstruction windows
 % param.N_elements : number of elements of the probe
 % param.width : pitch of individual transducers
 
 X_m = (1:param.N_elements)*(param.width) ;
 
 % M0 : image of origine (0,0) of the probe for each aquisition in order to
 % caracterize th initial rotation
 % X0 : x coordinate of the initial wavefront position
 % Z0 : x coordinate of the initial wavefront position (<0 for positive delays)
 
[theta,M0,X0,Z0] = EvalDelayLaw_shared(X_m,DelayLAWS,ActiveLIST,c); 


 Hf = figure(1);
 Ireconstruct = Retroprojection_shared( I , X_m , z_out , theta , M0 , Hf );

% subplot(211);plot(X_m*1e3,Ireconstruct(105,:)); subplot(212);plot(-20.5 + z_out*1e3,Ireconstruct(:,96))
 
% fwhm = FWHM(Ireconstruct(105,:),X_m*1e3)
% fwhm = FWHM(Ireconstruct(:,96),z_out*1e3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ho = figure(2) ;
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
cb = colorbar ;
title('Input phantom')
ylim([min(z_out*1e3) max(z_out*1e3)])
xlabel('x (mm)')
ylabel('z (mm)')
ylabel(cb,'light flux (a.u)')
colormap(parula)
set(findall(Ho,'-property','FontSize'),'FontSize',15)  



%rmpath('..\Scan routines')
%rmpath('functions')

