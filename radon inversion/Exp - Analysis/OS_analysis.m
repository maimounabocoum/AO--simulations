%% OF analysis %%
% addpath('..\shared functions folder');
%  MyImage = OS(Datas,ScanParam(:,1),ScanParam(:,2),...
%                 dFx,z,SampleRate*1e6,c,[X0 X1]*1e-3) ; 
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
X_m = (1:NbElemts)*(pitch*1e-3) ;
ElmtBorns   = [min(NbElemts,max(1,round(X0/pitch))),max(1,min(NbElemts,round(X1/pitch)))];
ElmtBorns   = sort(ElmtBorns) ; % in case X0 and X1 are mixed up
XMiddle  = mean(ElmtBorns)*(pitch*1e-3);

[ F_ct_kx , theta , decimation ] = MyImage.AddSinCos( MyImage.R ) ;
MyImage.F_R = MyImage.fourierz( F_ct_kx ) ;


    DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
    ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;   
    [theta,M0,~,~,C] = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;

    OriginIm = MyImage.iRadon( MyImage.F_R, X_m, XMiddle, MyImage.z , theta ,C , decimation , dFx);

    
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
%X_m = (1:128)*(0.11*1e-3) ;
[FTFx, theta , decimation ] = MyImage.AddSinCos(smooth(MyImage.R)) ;
    DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
    ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
       
    [theta,M0,~,~,C] = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;


FTF = MyImage.InverseFourierX( FTFx  , decimation , theta , C ) ;

OriginIm = sum(FTF,3) ;
figure;
imagesc(MyImage.x*1e3 + mean([X0 X1]),MyImage.z*1e3,real(OriginIm))
ylim([0 Prof]) 
xlim([X0,X1]) 
title('OS fourier')
% xlabel('x(mm)')
% ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u')
colormap(parula)
set(findall(gcf,'-property','FontSize'),'FontSize',15) 
