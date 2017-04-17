%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
field_init(0);

parameters;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);
% initial excitation field :
    Noc = 1; % number of optical cycles
    t_excitation = (0:1/param.fs:Noc*1.5/param.f0);
    excitation =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation))';    
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
% creating memory to save probe delay law
if param.Activated_FieldII == 1 
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);
end

Hf = figure(1);
for n_scan = 1:CurrentExperiement.Nscan

     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;
     % % option for screening : XY, Xt , XZt
     %   CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf)
     CurrentExperiement.MySimulationBox.ShowMaxField('XZ',Hf)
%      saveas(gcf,...
%      ['C:\Users\bocoum\Dropbox\PPT - prez\fieldII-images\OSmaxField_decimate',...
%      num2str(CurrentExperiement.ScanParam(n_scan,2)),'.png'],'png')
%      savefig(gcf,...
%      ['C:\Users\bocoum\Dropbox\PPT - prez\fieldII-images\OSmaxField_decimate',...
%      num2str(CurrentExperiement.ScanParam(n_scan,2)),'.fig'])
end
% CurrentExperiement.ShowAcquisitionLine();
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;