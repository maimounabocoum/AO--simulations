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

% Example : set active elements

%CurrentExperiement.MyProbe = CurrentExperiement.MyProbe.Set_ActiveList(1:300);

%% view active elements
CurrentExperiement.MyProbe.ShowProbe

%% set delay law
CurrentExperiement.MyProbe = ...
CurrentExperiement.MyProbe.Set_ActuatorDelayLaw('plane',-50*pi/180,param.c);

 
% view active elements
% CurrentExperiement.MyProbe.ShowDelay



