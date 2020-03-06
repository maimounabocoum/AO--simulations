%% OF analysis %%
 addpath('..\..\AO--commons\shared functions folder');
 

%%
% figure;
% imagesc(ScanParam(:,2),z*1e3,1e3*Datas)
% xlabel('order')
% ylabel('z(mm)')
% colormap(parula)
% cb = colorbar ;
% ylabel(cb, 'AC tension mV')
% set(findall(gcf,'-property','FontSize'),'FontSize',15) 

% [cx,cy,c] = improfile;
% figure;
% plot(sqrt((cx-cx(1)).^2 + (cy-cy(1)).^2),c)

%% reconstruction using iradon
% load : '2018-06-22'
c = 1540 ;
ElmtBorns   = [min(NbElemts,max(1,round(X0/pitch))),max(1,min(NbElemts,round(X1/pitch)))];
ElmtBorns   = sort(ElmtBorns) ; % in case X0 and X1 are mixed up
XMiddle  = mean(ElmtBorns)*(pitch*1e-3);

[ F_ct_kx , theta , decimation ] = MyImage.AddSinCos( MyImage.R ) ;
MyImage.F_R = MyImage.fourierz( F_ct_kx ) ;

    X_m = (1:param.N_elements)*(param.width) ;
    DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
    ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;   
    [theta,~,~,C] = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;

    
    OriginIm = MyImage.iRadon( MyImage.F_R, X_m, XMiddle, MyImage.z , theta ,C , decimation , dFx);

    
    Hresconstruct = figure;
   % set(Hresconstruct,'WindowStyle','docked');
    imagesc(X_m*1e3,MyImage.z*1e3,real(OriginIm));
    xlim([X0,X1])
    ylim([0 40])
    xlabel('x(mm)')
    ylabel('z(mm)')
    title('OS reconstruct')
    cb = colorbar;
    ylabel(cb,'a.u')
    colormap(parula)
    set(findall(Hresconstruct,'-property','FontSize'),'FontSize',15) 

     x = X_m*1e3;
     z = MyImage.z*1e3 ;
     save('C:\Users\mbocoum\Dropbox\self-written documents\acoustic-structured-illumination\images\datas\OPsimuIradon','x','z','OriginIm')

%% reconstruction using ifourier
c = 1540 ;
X_m = (1:NbElemts)*(pitch*1e-3) ;
[FTFx, theta , decimation ] = MyImage.AddSinCos(MyImage.R) ;
DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
       
[theta,~,~,C] = EvalDelayLawOS_shared( X_m-mean(X_m)  , DelayLAWS_  , ActiveLIST_ , c) ;


FTF = MyImage.InverseFourierX( FTFx  , decimation , theta , C ) ;

OriginIm = sum(FTF,3) ;
figure;
imagesc(MyImage.x*1e3 + mean([X0 X1]),MyImage.z*1e3,real(OriginIm))
ylim([0 30]) 
xlim([X0,X1]) 
title('OS fourier')
xlabel('x(mm)')
ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u')
colormap(parula)
set(findall(gcf,'-property','FontSize'),'FontSize',15) 

% save simulation for article
% x = MyImage.x*1e3 + mean([X0 X1]) ; 
% z = MyImage.z*1e3 ; 
% save('C:\Users\mbocoum\Dropbox\self-written documents\acoustic-structured-illumination\images\datas\OSsimuIFourier','FTF','x','z')
