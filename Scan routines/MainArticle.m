%% main script for article;

%% maimouna bocoum 04-01-2017
clearvars ;
addpath('..\Field_II')
addpath('..\radon inversion')
addpath('subscripts')
addpath('..\..\AO--commons\shared functions folder')
field_init(0);




 %% run this code portion to visualize the field temporal and/or spatial profile
 
 Hf = gcf;      % open a new figure
 n_scan =  1;  % index of the scan - look inside variable for corresponding parameters
 parameters; % script with simulation parameter (to edit befor running the simulation)
 CurrentExperiement = Experiment(param); 
 CurrentExperiement = CurrentExperiement.EvalPhantom();
 CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ;   % Initializes the Probe
 CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ;   % Calculate the Field Over input BOX

    % % option for screening : XY, Xt , XZt , Zt
     %CurrentExperiement.MySimulationBox.ShowMaxField('Xt', Hf);  
    % CurrentExperiement.MySimulationBox.ShowMaxField('Zt', Hf);  
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf);   
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZ', Hf);
   [Nx,Ny,Nz]       = SizeBox(CurrentExperiement.MySimulationBox); 
   x = CurrentExperiement.MySimulationBox.x ;
   z = CurrentExperiement.MySimulationBox.z ;
   
   t       = CurrentExperiement.MySimulationBox.time(:); % simulation time column vector
   myField = CurrentExperiement.GetCameraCorrelation(70e-6, 20e-6 ,n_scan);
   myField = reshape(myField ,[Ny,Nx,Nz]);     % resize the box to current screening
   myField = squeeze( myField(1,:,:) )' ; % remove single direction
                

              figure
     imagesc(x*1e3,z*1e3, myField );
     xlabel('x (mm)')
     ylabel('z (mm)')
     ylim([min(CurrentExperiement.MySimulationBox.z*1e3) max(CurrentExperiement.MySimulationBox.z*1e3)])
     title(['correlation with Reference, \tau_{exp,CCD} = ',num2str(20),'\mu s']) 
     cb = colorbar ;
     ylabel(cb,'a.u')
      drawnow
      
      % save parameters
     SubFolderName = generateSubFolderName('Q:\Data\simulations') ;
     FileName   = generateSaveName(SubFolderName ,'name','Profile','type',param.FOC_type,'Nbx',param.NbX,'Nbz',param.NbZ);
     save(FileName,'myField','x','z')  
             
%%


