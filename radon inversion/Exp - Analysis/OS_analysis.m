%% OF analysis %%
% addpath('..\shared functions folder');
% MyImage = OF(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
% MyImage = OP(x*1e-3,z*1e-3,Datas,SampleRate*1e6,c) ;
 MyImage = OS(Datas,ScanParam(:,1),ScanParam(:,2),...
                dFx,z,SampleRate*1e6,c,[X0 X1]*1e-3) ; 
% addpath('shared functions folder')

%%
figure;
imagesc(ScanParam(:,2),z*1e3,1e3*Datas)
xlabel('order')
ylabel('z(mm)')
colormap(parula)
cb = colorbar ;
ylabel(cb, 'AC tension mV')
set(findall(gcf,'-property','FontSize'),'FontSize',15) 

% [cx,cy,c] = improfile;
% figure;
% plot(sqrt((cx-cx(1)).^2 + (cy-cy(1)).^2),c)

%% reconstruction using iradon
X_m = (1:192)*(0.2*1e-3) ;


    MyImage.F_R = MyImage.fourierz(MyImage.R) ;   
    FILTER = MyImage.GetFILTER(2e-3);
    %MyImage.R   = MyImage.ifourierz((MyImage.F_R)) ;
    R   = MyImage.ifourierz((MyImage.F_R).*FILTER) ;
    
    [FTFx,~,decimation] = MyImage.AddSinCos(R) ;
   
    % resolution par iradon
    % FTF = MyImage.GetAngles(MyImage.R , decimation , theta ) ;
    % imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(FTF) );
    
    DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
    ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
       
    [theta,M0,~,~,C] = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;


    OriginIm = MyImage.Retroprojection( real(FTFx) , X_m, MyImage.z , theta , M0 , decimation , dFx);
%     OriginIm(1:400,:) = 0;
%     OriginIm(1000:1024,:) = 0;
    
    Hresconstruct = figure;
    set(Hresconstruct,'WindowStyle','docked');
    imagesc(X_m*1e3,MyImage.z*1e3,real(OriginIm));
    xlim([X0,X1])
    ylim([0 Prof])
    xlabel('x(mm)')
    ylabel('z(mm)')
    title('OS reconstruct')
    cb = colorbar;
    ylabel(cb,'a.u')
    colormap(parula)
    set(findall(Hresconstruct,'-property','FontSize'),'FontSize',15) 


%% reconstruction using ifourier
X_m = (1:192)*(0.2*1e-3) ;
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
% ylim([0 Prof]) 
% xlim([X0,X1]) 
OriginIm = OriginIm + FTF(:,:,nplot) ;
imagesc(MyImage.x*1e3,MyImage.z*1e3,real(OriginIm));   
T = unique(theta) ;
title(['theta = ',num2str(180*T(nplot)/pi)])
drawnow
pause(1)
end

figure;
imagesc(MyImage.x*1e3 + mean([X0 X1]),MyImage.z*1e3,real(OriginIm))
ylim([0 Prof]) 
xlim([X0,X1]) 
title('OS fourier')
xlabel('x(mm)')
ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u')
colormap(parula)
set(findall(gcf,'-property','FontSize'),'FontSize',15) 


%% view emission periods
figure;
imagesc(X_m*1e3,ScanParam(:,2),ActiveLIST')
colormap(parula)
