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
[angle, M0] = EvalDelayLawOS_shared( X_m , DelayLAWS, ActiveLIST, param.c);
 
for n_scan = 1:CurrentExperiement.Nscan
theta = angle(n_scan);
[Irad,Mcorner] = RotateTheta(X,Z,MyTansmission,theta);
d_offset = (Lprobe/2-M0(n_scan,1))*sin(theta) + (Z0-M0(n_scan,2))*cos(theta)-Z0;     
d_offset = d_offset/(z(2)-z(1)); % convert to pixels for sum 

Mask0 = interp1(CurrentExperiement.MyProbe.center(:,1,1)+ Lprobe/2,...
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
 Irad = Irad.*Mask0 ;
%Irad = Irad.*Mask ;  

Field_Profile(:,:,n_scan) = Mask0 ;

% correction matrice


imagesc(Irad)
AOSignal(:,n_scan) = interp1(1:size(Irad,1),trapz(x,Irad,2),...
                            (1:size(Irad,1))-d_offset,'linear',0) ;
drawnow

axis equal
end
 
AOSignal = AOSignal + 0*1e-3*rand(size(AOSignal)) ;
%%
figure;
imagesc(CurrentExperiement.ScanParam(:,2),z*1e3,AOSignal)
xlabel('order N_x')
ylabel('z(mm)')
% structure with appropriate axis and fourier transform methods
MyImage = OS(AOSignal,CurrentExperiement.ScanParam(:,1),...
             CurrentExperiement.ScanParam(:,2),param.df0x,...
             CurrentExperiement.MySimulationBox.z,...
             param.fs_aq,...
             param.c); 
          
MyImage.F_R = MyImage.fourierz( MyImage.R ) ; 
FILTER = MyImage.GetFILTER(1e-3);
MyImage.R   = MyImage.ifourierz(MyImage.F_R.*FILTER) ;

[FTFx, theta , decimation ] = MyImage.AddSinCos(MyImage.R) ;




%% resolution par iradon
% FTF = MyImage.GetAngles(MyImage.R , decimation , theta ) ;
DelayLAWS_  = MyImage.SqueezeRepeat(DelayLAWS) ;
ActiveLIST_ = MyImage.SqueezeRepeat( ActiveLIST ) ;

 c = 1540 ;
 
% %[theta,M0,X0,Z0]    = EvalDelayLawOS_shared( X_m , DelayLAWS_(:,1)  , ActiveLIST_(:,1) , c) ;
 [theta M0]    = EvalDelayLawOS_shared( X_m  , DelayLAWS_  , ActiveLIST_ , c) ;

 % Hf = figure;
 % X_m : interpolation vector for reconstruction
 % z :
 Ireconstruct = MyImage.Retroprojection( real(FTFx) , X_m , MyImage.z , theta , M0 , decimation , param.df0x);
%%
% FTFx : matrix with fourier composant : first cols = first decimation,
% vaying angle , second lines : second decimate, varying angle...

FTFxz = MyImage.fourierz(FTFx);
MyImage.ScatterFourier(FTFxz,decimation , theta);

FTF = MyImage.GetFourierX( FTFx  , decimation , theta ) ;

OriginIm = MyImage.ifourierx(FTF(:,:,1)) ;

figure('DefaultAxesFontSize',18); 
imagesc(MyImage.fx/MyImage.dfx,MyImage.z*1e3,abs(FTF(:,:,21)));

%imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(FTF));
%axis([-40 40 -100 100])
axis([-40 40 0 70])
title('reconstructed fourier along x')

figure('DefaultAxesFontSize',18);
imagesc(MyImage.x*1e3,MyImage.z*1e3,real(OriginIm));
xlabel('x(mm)')
ylabel('z(mm)')
xlim(param.Xrange*1000)
ylim(param.Zrange*1000)    
title('reconstructed object')


%% saving data to reconstruct folder
%test : fourier transform of original object
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 
%  figure('DefaultAxesFontSize',18);imagesc(x_phantom*1e3,z_phantom*1e3,MyTansmission)
 
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
 imagesc(MyImage.fx/MyImage.dfx,MyImage.z*1e3,abs(TinterpFFTx))
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
