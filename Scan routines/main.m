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
    excitation =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation))';
    %excitation =  exp(1i*(2*pi*param.f0*t_excitation -pi/2)).*hanning(length(t_excitation))'; % complexe field definition

    % figure;
    % plot(t_excitation*1e6,real(excitation),t_excitation*1e6,excitationRe)
    % xlabel('time in \mu s')
    % ylabel('a.u')
    % title('field excitation')
    
    % evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
%CurrentExperiement.ShowPhantom()    
 Hf = figure(1);
 tic
%h = waitbar(0,'Please wait...');
for n_scan = 1:CurrentExperiement.Nscan
%waitbar(n_scan/CurrentExperiement.Nscan)
     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf)
     
end
 


 toc
% % option for screening : XY, Xt , XZt
 CurrentExperiement.ShowAcquisitionLine(); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_end;