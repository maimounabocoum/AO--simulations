%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;
addpath('..\Field_II')
addpath('..\radon inversion')
addpath('subscripts')
addpath('..\..\AO--commons\shared functions folder')
field_init(0);
IsSaved = 0 ;

%%%%%%%%%%%% target folder to save simulated data %%%%%%%%%%
SimuPathFolder = 'Q:\datas\simulated datas';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% =========  %%%%%%%%%%%%%%%%%%%% Initialize Experiement  %%%%%%%%%%%%%%%%%%%%%%%
clearvars ; 

parameters; % script with simulation parameter (to edit befor running the simulation)
CurrentExperiement = Experiment(param); %initializes the experiement
  
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();

% use param.angles has an input to additionally show Radon transform
% CurrentExperiement.ShowPhantom(param.angles);
% CurrentExperiement.ShowPhantom();

% creating memory to save probe delay law
if param.Activated_FieldII == 1 
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);
end

 [Nx,Ny,Nz]    = SizeBox(CurrentExperiement.MySimulationBox);
 Field_Profile = zeros(Nz,Nx,CurrentExperiement.Nscan);
 
 %% run acquision loop over Nscan for tagged photon flux
 
 tic
 h = waitbar(0,'Please wait...');
 
 for n_scan = 1:CurrentExperiement.Nscan
 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ;   % Initializes the Probe
     CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ;   % Calculate the Field Over input BOX
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;   % Calculate current AO signal - Photorefractive
     % CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan,'Photorefractive','Holography') ; 
     
    % Retreive delay law for current scan for saving
    if strcmp(param.FOC_type,'OP') || strcmp(param.FOC_type,'OS')
     DelayLAWS( :  ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
    end
          
    waitbar(n_scan/CurrentExperiement.Nscan)
    
 end


 ActiveLIST = CurrentExperiement.BoolActiveList ; % list of active piezo-element for each n_scan (matrix), for saving
 
 close(h) 
 
 t_simulation = toc ;
 
 %% show acquisition loop results
 
 CurrentExperiement.ShowAcquisitionLine() ; % Show results of simulated acquistion
 % return the AO image on the simualtion grid for photorefractive detection
 % returns a plot for camera-based detection as resulted by camera
 % integration
 
 %% run this code portion to visualize the field temporal and/or spatial profile
 Hf = gcf;      % open a new figure
 n_scan =  1;  % index of the scan - look inside variable for corresponding parameters
 parameters; % script with simulation parameter (to edit befor running the simulation)
CurrentExperiement = Experiment(param); 
 CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ;   % Initializes the Probe
 CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ;   % Calculate the Field Over input BOX
     % % option for screening : XY, Xt , XZt
    % CurrentExperiement.MySimulationBox.ShowMaxField('Xt', Hf);  
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf);   
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZ', Hf);
    %
     %CurrentExperiement = CurrentExperiement.MyAO.BuildJM(param.f0,param.nuZ0,param.c,param.Bascule,CurrentExperiement.ScanParam);
     CurrentExperiement.ShowFieldCorrelation('XZ', Hf , 20e-6, 20e-6 ,n_scan); % ('XZ',Hf, startExposure, Exposure time,n_scan)
    % CurrentExperiement.MySimulationBox.ShowMaxField('YZ', Hf);
    
    
 %% run this code portion to visualize the field autocorrelation with reference beam (holography detection only)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% save data for reconstruction Iradon %% ONLY SAVING OP
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (IsSaved == 1)
     
     % saving folder name with todays date
     SubFolderName = generateSubFolderName(SimuPathFolder) ;
     FileName   = generateSaveName(SubFolderName ,'name','TaboltEffect_off','type',param.FOC_type);
     
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 y_phantom = CurrentExperiement.MySimulationBox.y ;
 z_phantom = CurrentExperiement.MySimulationBox.z ;
 [MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);
 
%-- saving parameters same as in experiment---------%

FreqSonde       = param.f0 ;
NbHemicycle     = 2*param.Noc ;
Foc             = param.focus ;
X0              = param.X0 ; 
X1              = param.X1 ;
pitch           = param.width ;
Prof            = max( z_phantom ) ;
ScanParam       = CurrentExperiement.ScanParam ;
NbElemts        = param.N_elements ;
SampleRate      = param.fs_aq ;
TypeOfSequence  = param.FOC_type ;
NbX             = param.NbX ;        % ondes JM
NbZ             = param.NbZ ;        % ondes JM
decimation      = param.decimation ; % ondes OS
dFx             = param.df0x;        % ondes OS
c               = param.c ;


%--------------------- saving datas -------------%
     switch param.FOC_type
         
         case 'OF'
 
 MyImage = OF(CurrentExperiement.ScanParam,CurrentExperiement.MySimulationBox.z,CurrentExperiement.AOSignal,param.fs_aq,param.c); 
 save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission','MyImage',...
               'DelayLAWS','ActiveLIST','MyImage','Field_Profile','param');  
 
         case 'OP'
 
%              AOSignal = CurrentExperiement.AOSignal ;
%              ScanParam = CurrentExperiement.ScanParam;

 MyImage = OP(CurrentExperiement.AOSignal,CurrentExperiement.ScanParam,...
              CurrentExperiement.MySimulationBox.z,param.fs_aq,param.c); 
            
 save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission','MyImage','R','zR',...
              'DelayLAWS','ActiveLIST','MyImage','Field_Profile','param');

         case 'OS'
  MyImage = OS(CurrentExperiement.AOSignal,CurrentExperiement.ScanParam(:,1),...
             CurrentExperiement.ScanParam(:,2),param.df0x,...
             CurrentExperiement.MySimulationBox.z,...
             param.fs_aq,...
             param.c,[param.X0 , param.X1]); 
     
    save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission',...
                  'DelayLAWS','ActiveLIST','MyImage','Field_Profile','param','Field')            
         case 'JM'  

    save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission',...
                  'DelayLAWS','ActiveLIST','param','CurrentExperiement')               
     end
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmpath('..\radon inversion')
% field_end;