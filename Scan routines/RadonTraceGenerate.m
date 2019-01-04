%%% generation of traces using radon transform for OP - OS parameters
% maimouna bocoum 24-10-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
addpath('..\radon inversion\shared functions folder')
field_init(0);

parameters;
IsSaved = 0 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
% I = CurrentExperiement.ShowPhantom();
[MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);

Lprobe = param.N_elements*(param.width + param.kerf) ;
Z0 = (max(CurrentExperiement.MySimulationBox.z)-min(CurrentExperiement.MySimulationBox.z))/2 ;


x = CurrentExperiement.MySimulationBox.x + Lprobe/2;
z = CurrentExperiement.MySimulationBox.z ;
[X,Z] = meshgrid(x,z);



% radon transform of the image :
AOSignal = zeros(length(z),CurrentExperiement.Nscan) ;
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);


 for n_scan = 1:CurrentExperiement.Nscan
theta = CurrentExperiement.ScanParam(n_scan,1);
CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan) ;
DelayLAWS( : ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
 end
 
X_m = (1:param.N_elements)*param.width; 
ActiveLIST = CurrentExperiement.BoolActiveList ;

% C : point of invariation by rotation of angle theta
[angles, M0 , ~ , ~ ,C] = EvalDelayLawOS_shared( X_m , DelayLAWS, ActiveLIST, param.c);
%% run scan 
figure ;
for n_scan = 1:CurrentExperiement.Nscan
% theta = angle(n_scan);
% C : center of rotation = [mean(X_m),0]
[Irad,MMcorner] = RotateTheta( X , Z , MyTansmission , angles(n_scan) , C(n_scan,:) );

u = [cos(theta) ; -sin(theta)] ;
v = [sin(theta) ; cos(theta)]  ;

Mask0 = interp1(CurrentExperiement.MyProbe.center(:,1,1)+ Lprobe/2 ,...
               double(CurrentExperiement.BoolActiveList(:,n_scan)),X);
           if n_scan ==1
           Mask = cos(0*X);    
           else
            
if mod(n_scan,2) == 0
Mask = cos(2*pi*param.df0x*CurrentExperiement.ScanParam(n_scan,2)*(X-Lprobe/2));
else
Mask = sin(2*pi*param.df0x*CurrentExperiement.ScanParam(n_scan,2)*(X-Lprobe/2));   
end
           end
% Irad = Irad.*Mask0 ;
Irad = Irad.*Mask ;  
Field_Profile(:,:,n_scan) = Mask0 ;

% correction matrice
imagesc(x*1e3,z*1e3,Irad)
xlabel('x(mm)')
ylabel('ct(mm)')
AOSignal(:,n_scan) = trapz(x,Irad,2) ;
% AOSignal(:,n_scan) = interp1(1:size(Irad,1),trapz(x,Irad,2),...
%                             (1:size(Irad,1))-d_offset,'linear',0) ;
drawnow
axis equal
end
 
AOSignal = AOSignal + 0*1e-3*rand(size(AOSignal)) ;
%%
% AOSignal = CurrentExperiement.AOSignal ;
% X_m = (1:param.N_elements)*param.width; 
% Lprobe = param.N_elements*(param.width + param.kerf) ;
% x = CurrentExperiement.MySimulationBox.x + Lprobe/2;
% z = CurrentExperiement.MySimulationBox.z ;

figure;
imagesc(CurrentExperiement.ScanParam(:,2),z*1e3,AOSignal)
xlabel('order N_x')
ylabel('z(mm)')
% structure with appropriate axis and fourier transform methods
MyImage = OS(AOSignal,CurrentExperiement.ScanParam(:,1),...
             CurrentExperiement.ScanParam(:,2),param.df0x,...
             CurrentExperiement.MySimulationBox.z,...
             param.fs_aq,...
             param.c,[min(X_m) , max(X_m)]);
   





%% resolution par iradon
[ F_ct_kx , theta , decimation ] = MyImage.AddSinCos( MyImage.R ) ;
MyImage.F_R = MyImage.fourierz( F_ct_kx ) ; 

% FTF = MyImage.GetAngles(MyImage.R , decimation , theta ) ;
DelayLAWS_  = MyImage.SqueezeRepeat( DelayLAWS  ) ;
ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;
[theta,M0,~,~,C]    = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , param.c) ;

 % Hf = figure;
 % X_m : interpolation vector for reconstruction
 % z :
 Ireconstruct = MyImage.iRadon( MyImage.F_R  , X_m , Lprobe/2, MyImage.z , theta , C , decimation , param.df0x);

 
figure
imagesc(X_m*1e3, MyImage.z*1e3,real(Ireconstruct))
xlim(param.Xrange*1000+ mean(X_m)*1000)
ylim(param.Zrange*1000) 
 
 %% iFourier
% FTFx : matrix with fourier composant : first cols = first decimation,
% vaying angle , second lines : second decimate, varying angle...

%FTFxz = MyImage.fourierz(FTFx);
% MyImage.ScatterFourier(FTFxz,decimation , theta);
[FTFx, theta , decimation ] = MyImage.AddSinCos(MyImage.R) ;
FTF = MyImage.InverseFourierX( FTFx  , decimation , theta , C ) ;
OriginIm = sum(FTF,3) ;


 
% OriginIm = 0 ;
% figure('DefaultAxesFontSize',18); 
% for nplot = 1:size(FTF,3)
% subplot(121)
% imagesc(MyImage.x*1e3,MyImage.z*1e3,real(FTF(:,:,nplot)));
% ylim(param.Zrange*1000) 
% xlim(param.Xrange*1000) 
% subplot(122)
% OriginIm = OriginIm + FTF(:,:,nplot) ;
% imagesc(MyImage.x*1e3,MyImage.z*1e3,real(OriginIm));   
% ylim(param.Zrange*1000) 
% xlim(param.Xrange*1000) 
% T = unique(theta) ;
% title(['theta = ',num2str(180*T(nplot)/pi)])
% drawnow
% pause(1)
% end

figure('DefaultAxesFontSize',18);
imagesc(MyImage.x*1e3 + mean(X_m)*1000,MyImage.z*1e3,real(OriginIm));
xlabel('x(mm)')
ylabel('z(mm)')
xlim(param.Xrange*1000+ mean(X_m)*1000)
ylim(param.Zrange*1000)    
title('reconstructed object')
colorbar


%% saving data to reconstruct folder
%test : fourier transform of original object
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 
%   figure('DefaultAxesFontSize',18);imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
 
 [X,Z] = meshgrid(MyImage.x,MyImage.z) ;
 [Xp,Zp] = meshgrid(x_phantom,z_phantom) ;
 
 Tinterp = interp2(Xp,Zp,MyTansmission,X,Z,'linear',0) ;
 TinterpFFT = MyImage.fourier( Tinterp );
 TinterpFFTx = MyImage.fourierx( Tinterp );
 figure('DefaultAxesFontSize',18);  
 imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(TinterpFFT))
 axis([-40 40 -100 100])
 xlabel('Fx/dfx')
 ylabel('Fz/dfz')
 title('object fourier transform')

  figure('DefaultAxesFontSize',18);  
 imagesc(MyImage.fx/MyImage.dfx,MyImage.z*1e3,angle(TinterpFFTx).*(abs(TinterpFFTx) > 0.01*max(abs(TinterpFFTx(:)))))
 axis([-40 40 0 70])
 xlabel('Fx/dfx')
 ylabel('z(mm)')
 title('object X fourier transform')

 
 if (IsSaved == 1)


     %--------------------- saving datas -------------
 save('..\radon inversion\saved images\SimulationTransmissionOS.mat','x_phantom','y_phantom','z_phantom','Field_Profile','MyTansmission','R','zR')
 save('..\radon inversion\saved images\SimulationOS.mat','MyImage','DelayLAWS','ActiveLIST','AOSignal')
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;
