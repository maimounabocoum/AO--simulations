%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
field_init(0);

parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);

% initial excitation field :

    Noc = 8; % number of optical cycles
    t_excitation = (0:1/param.fs:Noc*1.5/param.f0);
    excitation =  sin(2*pi*param.f0*t_excitation);
    excitation = excitation.*hanning(length(excitation))';

CurrentExperiement = CurrentExperiement.CalculateUSfield(excitation);



%CurrentExperiement = CurrentExperiement.StartExperiement();

% size(CurrentExperiement.MySimulationBox.Field)
% param.Nx*param.Ny*param.Nz
%
% option for screening : XY, Xt , XZt
CurrentExperiement.MySimulationBox.ShowMaxField('XZ'); % XZ : plan (x,z)
CurrentExperiement.ShowAcquisitionLine(); 
