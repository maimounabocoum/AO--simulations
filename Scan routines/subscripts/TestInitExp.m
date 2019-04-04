%% generation of traces using radon transform for OP - OS parameters
% maimouna bocoum 24-10-2017
clearvars ;

addpath('..\Field_II')
addpath('subscripts')
addpath('..\radon inversion')
addpath('..\radon inversion\shared functions folder')
field_init(0);


% edit parameters of the simulated experiment
parameters  ;


%% test initialization of an experiement

% creation of a class "Experiement". This class contains various
% sub-classes :
%  start inititialization of an experiment

CurrentExperiement = Experiment(param);
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
% I = CurrentExperiement.ShowPhantom();
[MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);

Lprobe = param.N_elements*(param.width + param.kerf) ;
Z0 = (max(CurrentExperiement.MySimulationBox.z)-min(CurrentExperiement.MySimulationBox.z))/2 ;
