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
    
clf;Hf = figure(1);
set(Hf,'NextPlot', 'replace');

 for n_scan = 1:CurrentExperiement.Nscan

     CurrentExperiement = CurrentExperiement.CalculateUSfield(excitation,n_scan);
     
    [Nx,Ny,Nz] = SizeBox(CurrentExperiement.MySimulationBox);
    [Field_max,Tmax] = max(CurrentExperiement.MySimulationBox.Field,[],1);
    Field_max = reshape(Field_max,[Ny,Nx,Nz]);
       
    CurrentExperiement.MyProbe.ShowProbe();

           close all
           imagesc(CurrentExperiement.MySimulationBox.x*1e3,CurrentExperiement.MySimulationBox.z*1e3,...
                    squeeze(Field_max(1,:,:))');
   
            shading interp
            xlabel('x (mm)')
            ylabel('z (mm)')
            title(['Maximum Field in plane Y = ',num2str(CurrentExperiement.MySimulationBox.y(1)*1e3),'mm'])
            colorbar
            drawnow
 end
 
% % option for screening : XY, Xt , XZt
% 
% CurrentExperiement.ShowAcquisitionLine(); 
