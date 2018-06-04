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
% MyImage = OP(data(:,:,1),X*pi/180,Y*1e-3,Param.SamplingRate*1e6,c); 

%% reconstruction using iradon
X_m = (1:192)*(0.2*1e-3) ;


    MyImage.F_R = MyImage.fourierz(MyImage.R) ;   
    FILTER = MyImage.GetFILTER(1e-3);
    %MyImage.R   = MyImage.ifourierz((MyImage.F_R)) ;
    R   = MyImage.ifourierz((MyImage.F_R).*FILTER) ;
    
    [FTFx,~,decimation] = MyImage.AddSinCos(R) ;
   
    % resolution par iradon
    % FTF = MyImage.GetAngles(MyImage.R , decimation , theta ) ;
    % imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(FTF) );
    
    DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
    ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
       
    [theta,M0,~,~,C] = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;


    OriginIm = MyImage.Retroprojection( real(FTFx) , X_m, MyImage.z , theta , M0 , decimation ,MyImage.dfx);
%     OriginIm(1:400,:) = 0;
%     OriginIm(1000:1024,:) = 0;
    
    Hresconstruct = figure;
    set(Hresconstruct,'WindowStyle','docked');
    imagesc(X_m*1e3,MyImage.z*1e3,real(OriginIm));
    xlim(MyImage.Lx*1e3+20)
    ylim(MyImage.Lz*1e3)
    xlabel('x(mm)')
    ylabel('z(mm)')
    title('OS reconstruct')
    cb = colorbar;
    ylabel(cb,'a.u')
    colormap(parula)
    set(findall(Hresconstruct,'-property','FontSize'),'FontSize',15) 


%% inversion through inverse fourier transform:

%% reconstruction using ifourier
X_m = (1:192)*(0.2*1e-3) - mean(X_m );
[FTFx,~,decimation] = MyImage.AddSinCos(MyImage.R) ;
DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
[theta,M0,~,~,C] = EvalDelayLawOS_shared( X_m , DelayLAWS_  , ActiveLIST_ , c) ;

FTF = MyImage.InverseFourierX( FTFx  , decimation , theta , C ) ;



    
OriginIm = 0 ;
figure; 
for nplot = 1:size(FTF,3)
subplot(121)
imagesc(MyImage.x*1e3,MyImage.z*1e3,real(FTF(:,:,nplot)));
ylim(MyImage.Lz*1e3)
% xlim([X0,X1]) 
OriginIm = OriginIm + FTF(:,:,nplot) ;
imagesc(MyImage.x*1e3,MyImage.z*1e3,real(OriginIm));   
T = unique(theta) ;
title(['theta = ',num2str(180*T(nplot)/pi)])
drawnow
pause(1)
end

figure;
imagesc(MyImage.x*1e3+20,MyImage.z*1e3,real(OriginIm))
%     xlim(MyImage.Lx*1e3+20)
ylim(MyImage.Lz*1e3)
title('OS fourier')
xlabel('x(mm)')
ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u')
colormap(parula)
set(findall(gcf,'-property','FontSize'),'FontSize',15) 


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


