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
  load('saved images\SimulationOS_field.mat');
  load('saved images\SimulationTransmissionOS.mat');
  load('saved images\Field_Profile.mat');
 %% test reconstruction
%  
%  Iref = sum(Field_Profile,3);
%  figure; imagesc(Iref);
FTF = MyImage.AddSinCosRef(MyImage.rawData, Field_Profile) ;
figure; imagesc(abs(FTF)) 
 
 %% fourier reconstruction
MyImage.F_R = MyImage.fourierz( MyImage.R ) ;    
[MyImage.F_R, MyImage.theta,MyImage.decimation ] = MyImage.AddSinCos(MyImage.F_R) ;
FTF = MyImage.GetFourier(MyImage.F_R,MyImage.decimation ) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OriginIm = MyImage.ifourier(FTF) ;

figure('DefaultAxesFontSize',18); 
imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(FTF) );
axis([-40 40 -100 100])
title('reconstructed fourier transform')

figure; imagesc( abs(OriginIm) );
title('reconstructed object')
%figure;imagesc(CurrentExperiement.BoolActiveList)

%% saving data to reconstruct folder
%test : fourier transform of original object
 [X,Z] = meshgrid(MyImage.x,MyImage.z) ;
 [Xp,Zp] = meshgrid(x_phantom,z_phantom) ;
 
 Tinterp = interp2(Xp,Zp,MyTansmission,X,Z,'linear',0) ;
 TinterpFFT = MyImage.fourier( Tinterp );
 
 figure('DefaultAxesFontSize',18);  
 imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(TinterpFFT))
 axis([-40 40 -100 100])
 xlabel('Fx/dfx')
 ylabel('Fz/dfz')
 title('object fourier transform')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting the final results and its fourier transform
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('DefaultAxesFontSize',18);  
imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
colorbar
title('simulation input phantom')
xlabel('x (mm)')
ylabel('y (mm)')
drawnow   


