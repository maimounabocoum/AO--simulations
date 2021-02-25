%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;
addpath('..\Field_II')
addpath('..\radon inversion')
addpath('subscripts')
addpath('..\..\AO--commons\shared functions folder')
field_init(0);





%% =========  %%%%%%%%%%%%%%%%%%%% Initialize Experiement  %%%%%%%%%%%%%%%%%%%%%%%
clearvars ; 
IsSaved = 1 ;
%%%%%%%%%%%% target folder to save simulated data %%%%%%%%%%
SimuPathFolder = 'D:\Data\simulations';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



parameters; % script with simulation parameter (to edit befor running the simulation)
CurrentExperiement = Experiment(param); %initializes the experiement
  
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();

% use param.angles has an input to additionally show Radon transform
% CurrentExperiement.ShowPhantom(param.angles);
 Iphantom = CurrentExperiement.ShowPhantom();

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
 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ;       % Initializes the Probe
     CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ;       % Calculate the Field Over input BOX
     %CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;      % Calculate current AO signal - Photorefractive
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan,'camera') ;      % detection choice : photorafractive / holography
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan,'photodiode') ;  % detection choice : photorafractive / holography
     
    % Retreive delay law for current scan for saving
    if strcmp(param.FOC_type,'OP') || strcmp(param.FOC_type,'OS')
     DelayLAWS( :  ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
    end
          
    waitbar(n_scan/CurrentExperiement.Nscan)
    
 end


 ActiveLIST = CurrentExperiement.BoolActiveList ; % list of active piezo-element for each n_scan (matrix), for saving
 
 close(h) 
 
 t_simulation = toc ;
 
 %% show acquisition loop results CCD
 
 n_scan = 3 ; % choose which scan to screen out
 [Nx,Ny,Nz] = CurrentExperiement.MySimulationBox.SizeBox();
 MySignal  = CurrentExperiement.AOSignal_CCD(:,n_scan);
 MySignal = reshape(MySignal,[Ny,Nx,Nz]);
 MySignal = squeeze( MySignal(1,:,:) )' ; % remove single direction
 figure; imagesc(MySignal);
 
 %% show acquisition loop results photodiode
 
 CurrentExperiement.ShowAcquisitionLine('photodiode') ; % Show results of simulated acquistion
 CurrentExperiement.ShowAcquisitionLine('camera') ; % Show results of simulated acquistion

 %% view fourier matrix for specific comonant
 n_scan = 10;
  CurrentExperiement.ShowActualFFT(n_scan) ; % Show results of simulated acquistion


 %% reconstruction of JM image:
 
 CurrentExperiement.ShowJMreconstruction_photodiode();%inpout phase

 %%
 
  CurrentExperiement.ShowJMreconstruction_camera();
  
  %% matrix reconstruction
  % Iphantom = CurrentExperiement.ShowPhantom();
  [Nx,Ny,Nz] = CurrentExperiement.MySimulationBox.SizeBox();
  G = JM( Nx , Nz , Nx*param.nuX0 , Nz*param.nuZ0 );
  
  [ Y , yNscan , yScanParam ] = CurrentExperiement.GetYvector('Real'); 
  %M = CurrentExperiement.GetMmatrix('Real'); % 'Real' , 'Real-4phase','Fourier','Fourier-4phase'
  
  [G,Iout] = G.BuiltRealMatrix( M , CurrentExperiement.ScanParam , -5:5 , 1:10 , [0,0.25,0.5,0.75]);
  
Y = Y(Iout);  
[U,S,V] = svd(G.M);
 Splus = 0*S;
 seuil = S(240,240)/4;
 Splus(S>seuil) = 1./S(S>seuil);
 Splus = Splus';
 Minv = V*Splus*U';
  Ireconst = Minv*Y ;
  
figure;
subplot(121)
imagesc(Iphantom)
subplot(122)
imagesc( squeeze(reshape( Ireconst ,[Ny,Nx,Nz]))' )

 
  %%
  
%   for i=1:size(M,1)
%      Fplane = CurrentExperiement.vector2FFTplane( M(i,:) , G );
%      imagesc(G.fx/(G.dfx),G.fz/(G.dfz),abs(Fplane))
%      drawnow
%   end
  
  % M : 'Fourier-4phase'
  % each line is the fourier transform of g-illumination pattern
 

  
  % reconstruct image from Matrix:
  figure(1)
  subplot(121)
      FplaneIm = CurrentExperiement.GetYvector('Real-4phase');
      Fplane = CurrentExperiement.vector2FFTplane(FplaneIm, G );
      imagesc(G.fx/(G.dfx),G.fz/(G.dfz),abs(Fplane))
      colorbar
      title('fetch in experiment')
  subplot(122)
      FplaneImFFT = CurrentExperiement.GetYvector('Fourier-4phase');
      Fplane = CurrentExperiement.vector2FFTplane( M*FplaneImFFT , G );
      imagesc(G.fx/(G.dfx),G.fz/(G.dfz),abs(Fplane))
      colorbar
      title('calculated using matrix approach')
      
%       figure;
%       plot( abs(FplaneIm) )
%       hold on
%       plot( abs(M*FplaneImFFT) )
   
%%

  CurrentExperiement.ShowMatrixreconstruction_camera()

  
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
    % myField = CurrentExperiement.GetCameraTagged(20e-6,20e-6,n_scan);
    % CurrentExperiement.MyAO = CurrentExperiement.MyAO.AOsequenceGenerate(param,CurrentExperiement.ScanParam);
   %FieldCorr =  CurrentExperiement.GetCameraCorrelation(20e-6, 20e-6 ,n_scan);
      CurrentExperiement.ShowFieldCorrelation('XZ', Hf , 20e-6, 20e-6 ,n_scan); % ('XZ',Hf, startExposure, Exposure time,n_scan)
    % CurrentExperiement.MySimulationBox.ShowMaxField('YZ', Hf);
%  CurrentExperiement.MyProbe.ShowProbe()   
% figure;imagesc(CurrentExperiement.MySimulationBox.Field)    
 %% run this code portion to visualize the field autocorrelation with reference beam (holography detection only)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% save data for reconstruction Iradon %% ONLY SAVING OP
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (IsSaved == 1)
     
     % saving folder name with todays date
     SubFolderName = generateSubFolderName(SimuPathFolder) ;
     FileName   = generateSaveName(SubFolderName ,'name','Simu_longwindows_NbZ=10_NbX=10','type',param.FOC_type);
     
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
phase           = param.phase  ;     % ondes JM
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
                  'param','CurrentExperiement')               
     end
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmpath('..\radon inversion')
% field_end;