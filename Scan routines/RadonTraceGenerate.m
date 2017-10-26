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
I = CurrentExperiement.ShowPhantom();

Lprobe = param.N_elements*(param.width + param.kerf) ;
Z0 = (max(CurrentExperiement.MySimulationBox.z)-min(CurrentExperiement.MySimulationBox.z))/2 ;


x = CurrentExperiement.MySimulationBox.x + Lprobe/2;
z = CurrentExperiement.MySimulationBox.z ;
[X,Z] = meshgrid(x,z);



% radon transform of the image :
AOSignal = zeros(length(z),CurrentExperiement.Nscan) ;
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);


 for n_scan = 1:CurrentExperiement.Nscan
theta = CurrentExperiement.ScanParam(n_scan);
CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan) ;
DelayLAWS( CurrentExperiement.MyProbe.ActiveList ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
 end
 
X_m = (0:191)*0.2*1e-3 ; 
ActiveLIST = CurrentExperiement.BoolActiveList ;
[angle, M0] = EvalDelayLawOS_shared( X_m , DelayLAWS , ActiveLIST , param.c);
 

for n_scan = 1:CurrentExperiement.Nscan
theta = angle(n_scan);
[Irad,Mcorner] = RotateTheta(X,Z,I,theta);
%d_offset = ( R*tan(CurrentExperiement.ScanParam(n_scan,1)/2) )*tan(CurrentExperiement.ScanParam(n_scan,1));
d_offset = (Lprobe/2-M0(n_scan,1))*sin(theta) + (Z0-M0(n_scan,2))*cos(theta)-Z0;     
dtest(n_scan) = d_offset;
% d_offset = abs(tan(theta)*Lprobe/2 - R - tan(theta)*(Lprobe-sign(theta)*Lprobe)/2)/sqrt(1+tan(theta)^2) - R;     
d_offset = d_offset/(z(2)-z(1)); % convert to pixels for sum 
Mask = interp1(CurrentExperiement.MyProbe.center(:,1,1),double(CurrentExperiement.BoolActiveList(:,n_scan)),X);

%Irad = Irad.*Mask ;
%  imagesc(Irad)
AOSignal(:,n_scan) = interp1(1:size(Irad,1),sum(Irad,2),(1:size(Irad,1))-d_offset,'linear',0) ;
% drawnow
% pause(0.2)
% axis equal
 end

figure;
imagesc(CurrentExperiement.ScanParam*180/pi,z*1e3,AOSignal)
MyImage = OP(AOSignal,CurrentExperiement.ScanParam,CurrentExperiement.MySimulationBox.z,param.fs_aq,param.c); 
 

%% saving data to reconstruct folder

 if (IsSaved == 1)
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 [MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);
 
     %--------------------- saving datas -------------
 save('..\radon inversion\saved images\SimulationTransmissionOS.mat','x_phantom','y_phantom','z_phantom','MyTansmission','R','zR')
     if param.Activated_FieldII == 1
     save('..\radon inversion\saved images\SimulationOS.mat','MyImage','DelayLAWS','ActiveLIST')
     else
     save('..\radon inversion\saved images\SimulationOS.mat','MyImage')
     end
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;
