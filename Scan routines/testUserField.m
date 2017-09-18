%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
field_init(0);

parameters;
IsSaved = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);

% set active profile 
CurrentExperiement.MyProbe = ...
CurrentExperiement.MyProbe.Set_ActiveList(40:100);

% set delay law
CurrentExperiement.MyProbe = ...
CurrentExperiement.MyProbe.Set_ActuatorDelayLaw('plane',-50*pi/180,param.c);

 
% MyField = ExcitationField( CurrentExperiement.MyProbe , param.f0 , param.fs , param.Noc );
% MyField = MyField.Propagate((1:100)*1e-3,1540);




