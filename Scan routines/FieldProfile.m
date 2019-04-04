%%% generation of field spatial profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
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

% initial excitation field :

    t_excitation = (0:1/param.fs:param.Noc*1.5/param.f0);
    excitation   =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation)).^2';
    
%     excitation_env = hilbert(excitation);
%     excitation_env= abs(excitation_env);
% 
%     figure;
%     plot(t_excitation*1e6,excitation)
%     hold on 
%     plot(t_excitation*1e6,excitation_env,'color','red')
%     xlabel('time in \mu s')
%     ylabel('a.u')
%     title('field excitation')
    

    
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();
CurrentExperiement.ShowPhantom();
%use param.angles has an input to additionally show Radon transform


% creating memory to save probe delay law
if param.Activated_FieldII == 1 
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);
end

 tic
 Hf = gcf;
 h = waitbar(0,'Please wait...');

 [Nx,Ny,Nz]   = SizeBox(CurrentExperiement.MySimulationBox);
 Field_Profile = zeros(Nz,Nx,CurrentExperiement.Nscan);
 
 for n_scan = 1%:CurrentExperiement.Nscan
 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan);
     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     % % option for screening : XY, Xt , XZt

     CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf)    
    %  CurrentExperiement.MySimulationBox.ShowMaxField('XZ', Hf)

            [Field_max,Tmax] = max(CurrentExperiement.MySimulationBox.Field,[],1);
            % max(obj.Field,[],1) : returns for each colulm
            % the maximum field pressure.
            Field_Profile(:,:,n_scan) = squeeze( reshape(Field_max,[Ny,Nx,Nz]) )';

            
    % retreive delay law for current scan
    if strcmp(param.FOC_type,'OP') || strcmp(param.FOC_type,'OS')
     DelayLAWS( :  ,n_scan) = ...
                CurrentExperiement.MyProbe.DelayLaw ;
    end
          
    waitbar(n_scan/CurrentExperiement.Nscan)
    
 end
 
 %% screening loop
%  for n_scan = 1:CurrentExperiement.Nscan
%     imagesc(Field_Profile(:,:,n_scan))
%     caxis([min(min(Field_Profile(:,:,n_scan)))...
%            max(max(Field_Profile(:,:,n_scan))) ])
%        drawnow
%        pause(0.5)
%  end
%  
 