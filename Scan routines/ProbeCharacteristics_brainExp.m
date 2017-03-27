%%% Study of probe characterisics %%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
field_init(0);

parameters;


    Noc = 1; % number of optical cycles
    t_excitation = (0:1/param.fs:Noc*1.5/param.f0);
    excitation =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation))';

        Hf = figure(1);

% initializetion of parametetric study
focus = 8/1000;


Field2D = [];
for loop = 1:length(focus)
param.focus = focus(loop);
% initial excitation field :

     CurrentExperiement = Experiment(param);
     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,64);
     %   CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;
     % % option for screening : XY, Xt , XZt
     
            [Field_max,Tmax] = max(CurrentExperiement.MySimulationBox.Field,[],1);
            % max(obj.Field,[],1) : returns for each colulm
            % the maximum field pressure.
            Field_max = reshape(Field_max,[param.Ny,param.Nx,param.Nz]);
            Tmax = reshape(Tmax,[param.Ny,param.Nx,param.Nz]);
           
            Field2D = [Field2D,squeeze(Field_max(1,:,:))'];
            
            % extraction of z profil :
            
            subplot(1,length(focus),loop)
            imagesc(CurrentExperiement.MySimulationBox.x*1e3,CurrentExperiement.MySimulationBox.z*1e3,Field2D);
            xlabel('x (mm)')
            ylabel('z (mm)')
            colorbar
            drawnow
            
            
            % resolution along Z:
             figure;
             plot(CurrentExperiement.MySimulationBox.z*1e3,Field2D(:,74))
             fhwm( CurrentExperiement.MySimulationBox.z*1e3 ,Field2D(:,74) )
%             hold on 
%             plot(CurrentExperiement.MySimulationBox.x*1e3,Field2D(49,:))
%             fhwm(CurrentExperiement.MySimulationBox.x*1e3,Field2D(49,:))

end

% Xim = 1:size(Field2D,2);
% imagesc(Xim,CurrentExperiement.MySimulationBox.z*1e3,Field2D)
% xlabel('focus = (5:0.5:10)/1000;')
% ylabel('z (mm)')

field_end;












