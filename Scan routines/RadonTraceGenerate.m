%%% generation of traces using radon transform for OP - OS parameters
% maimouna bocoum 24-10-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
addpath('..\radon inversion\shared functions folder')
field_init(0);

parameters;
IsSaved = 1 ;

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
DelayLAWS( CurrentExperiement.MyProbe.ActiveList ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
 end
 
X_m = (0:191)*0.2*1e-3 ; 
ActiveLIST = CurrentExperiement.BoolActiveList ;
[angle, M0] = EvalDelayLawOS_shared( X_m , DelayLAWS , ActiveLIST , param.c);
 
for n_scan = 1:CurrentExperiement.Nscan
theta = angle(n_scan);
[Irad,Mcorner] = RotateTheta(X,Z,MyTansmission,theta);
d_offset = (Lprobe/2-M0(n_scan,1))*sin(theta) + (Z0-M0(n_scan,2))*cos(theta)-Z0;     
d_offset = d_offset/(z(2)-z(1)); % convert to pixels for sum 
Mask0 = interp1(CurrentExperiement.MyProbe.center(:,1,1)+ Lprobe/2,...
               double(CurrentExperiement.BoolActiveList(:,n_scan)),X);
if mod(n_scan,2) == 0
Mask = sin(2*pi*param.df0x*CurrentExperiement.ScanParam(n_scan,2)*X);
else
Mask = cos(2*pi*param.df0x*CurrentExperiement.ScanParam(n_scan,2)*X);   
end

Irad = Irad.*Mask0 ;
%Irad = Irad.*Mask ;
% plot(x,Mask0)
% hold on
% plot(x,Mask)
% hold off
imagesc(Irad)
AOSignal(:,n_scan) = interp1(1:size(Irad,1),trapz(x,Irad,2),...
                            (1:size(Irad,1))-d_offset,'linear',0) ;
drawnow


axis equal
end
 
figure;
imagesc(CurrentExperiement.ScanParam(:,1)*180/pi,z*1e3,AOSignal)

MyImage = OS(AOSignal,CurrentExperiement.ScanParam(:,1),...
             CurrentExperiement.ScanParam(:,2),param.df0x,...
             CurrentExperiement.MySimulationBox.z,...
             param.fs_aq,...
             param.c); 
          
MyImage.F_R = MyImage.fourierz( MyImage.R ) ;    
[MyImage.F_R, MyImage.theta,MyImage.decimation ] = MyImage.AddSinCos(MyImage.F_R) ;
FTF = MyImage.GetFourier(MyImage.F_R,MyImage.decimation ) ;
 
figure; imagesc( abs(FTF) );
 
 % test : fourier transform of original object
%  [X,Z] = meshgrid(MyImage.x,MyImage.z) ;
%  [Xp,Zp] = meshgrid(x_phantom,z_phantom) ;
%  
%  Tinterp = interp2(Xp,Zp,MyTansmission,X,Z,'linear',0) ;
%  TinterpFFT = MyImage.fourier( Tinterp );
%  figure; imagesc(MyImage.fx/MyImage.dfx,MyImage.fz/MyImage.dfz,abs(TinterpFFT))

%% inverse fourier transform on calculated datas :
% [Angles,~,Iangles] = unique(obj.theta) ;
% for iangle=1:length(Angles)
% Iextract = find(Iangles == iangle);

OriginIm = MyImage.ifourier(FTF) ;
figure; imagesc( abs(fftshift(OriginIm,2)) );
figure; imagesc( abs(OriginIm) );


%% saving data to reconstruct folder

 if (IsSaved == 1)
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;

     %--------------------- saving datas -------------
 save('..\radon inversion\saved images\SimulationTransmissionOS.mat','x_phantom','y_phantom','z_phantom','MyTansmission','R','zR')
 save('..\radon inversion\saved images\SimulationOS.mat','MyImage','DelayLAWS','ActiveLIST','AOSignal')
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;
