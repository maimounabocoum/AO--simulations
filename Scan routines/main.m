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
%     figure;
%     plot(t_excitation*1e6,real(excitation))
%     xlabel('time in \mu s')
%     ylabel('a.u')
%     title('field excitation')
    
    % evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
%use param.angles has an input to additionally show Radon transform


% creating memory to save probe delay law
if param.Activated_FieldII == 1 
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);
end

 Hf = figure(1);
 tic
h = waitbar(0,'Please wait...');
for n_scan = 1:CurrentExperiement.Nscan
 waitbar(n_scan/CurrentExperiement.Nscan)

     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;
     % % option for screening : XY, Xt , XZt

    %CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf)
     
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZ',Hf)
   
    % retreive delay law for cuurent scan
    if param.Activated_FieldII == 1
     DelayLAWS(:,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
    end
    
       
      
end
 
 toc
 CurrentExperiement.ShowAcquisitionLine();
 
% [Nx,Ny,Nz] = CurrentExperiement.MySimulationBox.SizeBox();
% Transmission = squeeze( reshape(CurrentExperiement.DiffuseLightTransmission',[Ny,Nx,Nz]) );
% plot(CurrentExperiement.MySimulationBox.z*1e3,Transmission(75,:))
%hold on
%plot(CurrentExperiement.MySimulationBox.z*1e3,CurrentExperiement.AOSignal(:,64)/max(CurrentExperiement.AOSignal(:,64)))
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% save data for reconstruction Iradon %% ONLY SAVING OP
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if param.FOC_type == 'OP'
 MyImage = OP(CurrentExperiement.AOSignal,CurrentExperiement.ScanParam,CurrentExperiement.MySimulationBox.z,param.fs_aq,param.c); 
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 [MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);
 
     %--------------------- saving datas -------------
 save('..\radon inversion\saved images\SimulationTransmission.mat','x_phantom','y_phantom','z_phantom','MyTansmission','R','zR')
     if param.Activated_FieldII == 1
     save('..\radon inversion\saved images\Simulation.mat','MyImage','DelayLAWS')
     else
     save('..\radon inversion\saved images\Simulation.mat','MyImage')
     end
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmpath('..\radon inversion')
field_end;