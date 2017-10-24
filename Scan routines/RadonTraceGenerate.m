%%% generation of traces using radon transform for OP - OS parameters
% maimouna bocoum 24-10-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
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
[X,Z] = meshgrid(CurrentExperiement.MySimulationBox.x,CurrentExperiement.MySimulationBox.z);



% radon transform of the image :
AOSignal = zeros(length(CurrentExperiement.MySimulationBox.z),CurrentExperiement.Nscan) ;

 for n_scan = 1:CurrentExperiement.Nscan

Irad = imrotate(I,CurrentExperiement.ScanParam(n_scan,1)*180/pi,'loose','crop');

Mask = interp1(CurrentExperiement.MyProbe.center(:,1,1),double(CurrentExperiement.BoolActiveList(:,n_scan)),X);

Irad = Irad.*Mask ;
imagesc(CurrentExperiement.MySimulationBox.x,CurrentExperiement.MySimulationBox.z,Irad)
AOSignal(:,n_scan) = sum(Irad,2) ;

 end

 figure;
imagesc(AOSignal)


MyImage = OP(AOSignal,ScanParam,MySimulationBox.z,fs_aq,c); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;
