%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
addpath('subscripts')
addpath('..\..\AO--commons\shared functions folder')
field_init(0);

parameters;
IsSaved = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);

% initial excitation field :

%     excitation_env = hilbert(excitation);
%     excitation_env= abs(excitation_env);
% 
%     figure;
%     plot(t_excitation*1e6,excitation)
%     hold on 
%     plot(t_excitation*1e6,excitation,'color','red')
%     xlabel('time in \mu s')
%     ylabel('a.u')
%     title('field excitation')
    
    
% evaluate Phantom on simulation Box :
CurrentExperiement = CurrentExperiement.EvalPhantom();

%CurrentExperiement.ShowPhantom();
%use param.angles has an input to additionally show Radon transform


% creating memory to save probe delay law
if param.Activated_FieldII == 1 
DelayLAWS = zeros(param.N_elements,CurrentExperiement.Nscan);
end

 [Nx,Ny,Nz]    = SizeBox(CurrentExperiement.MySimulationBox);
 Field_Profile = zeros(Nz,Nx,CurrentExperiement.Nscan);
 
 %% run acquision loop over Nscan
 
% figure(1);plot(real(CurrentExperiement.MyAO.Event)); 

 tic
 Hf = gcf;
 h = waitbar(0,'Please wait...');

 

 
 for n_scan = 1:CurrentExperiement.Nscan
 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ; % Initializes the Probe
     CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ; % Calculate the Field Over input BOX
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ; % Calculate current AO signal - Photorefractive
     %  CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan,'Photorefractive','Holography') ; 
     % Calculate current AO signal - Holographie

     

    % % option for screening : XY, Xt , XZt
    % CurrentExperiement.MySimulationBox.ShowMaxField('Xt', Hf);  
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf);   
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZ', Hf);
    % CurrentExperiement.ShowFieldCorrelation('XZ', Hf , 23.5e-6 , 1);
    % CurrentExperiement.MySimulationBox.ShowMaxField('YZ', Hf);
    
    % field profile
    [Field_max,Tmax] = max(CurrentExperiement.MySimulationBox.Field,[],1);
    % max(obj.Field,[],1) : returns for each colulm
    % the maximum field pressure.
    % Field_Profile(:,:,n_scan) = squeeze( reshape(Field_max,[Ny,Nx,Nz]) )';
   
    % retreive delay law for cuurent scan
    if strcmp(param.FOC_type,'OP') || strcmp(param.FOC_type,'OS')
     DelayLAWS( :  ,n_scan) = CurrentExperiement.MyProbe.DelayLaw ;
    end
          
    waitbar(n_scan/CurrentExperiement.Nscan)
    
 end


 ActiveLIST = CurrentExperiement.BoolActiveList ;
 
 close(h) 
 
 t_simulation = toc
 
 %% show acquisition loop results
 
 
 CurrentExperiement.ShowAcquisitionLine();
 
%  figure
%  imagesc(CurrentExperiement.ScanParam*1e3+20,...
%           CurrentExperiement.MySimulationBox.z*1e3,...
%           CurrentExperiement.AOSignal)
%       
%   x = CurrentExperiement.ScanParam*1e3+20 ;
%   z = CurrentExperiement.MySimulationBox.z*1e3 ;
%   I = CurrentExperiement.AOSignal ;
%   save('C:\Users\mbocoum\Dropbox\self-written documents\acoustic-structured-illumination\images\datas\SimuOF','x','z','I')

% [Nx,Ny,Nz] = CurrentExperiement.MySimulationBox.SizeBox();
% Transmission = squeeze( reshape(CurrentExperiement.DiffuseLightTransmission',[Ny,Nx,Nz]) );
% plot(CurrentExperiement.MySimulationBox.z*1e3,Transmission(75,:))
% hold on
% plot(CurrentExperiement.MySimulationBox.z*1e3,CurrentExperiement.AOSignal(:,64)/max(CurrentExperiement.AOSignal(:,64)))
 
% set(findall(gcf,'-property','FontSize'),'FontSize',15) 
% [cx,cy,c] = improfile;
% figure;
% plot(cx(1) + sqrt((cx-cx(1)).^2 + (cy-cy(1)).^2),c/max(c))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% save data for reconstruction Iradon %% ONLY SAVING OP
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if (IsSaved == 1)
     
     % saving folder name with todays date
     SubFolderName = generateSubFolderName('..\radon inversion\saved images') ;
     FileName   = generateSaveName(SubFolderName ,'name','500micronsHole_Foc19mm','type',param.FOC_type);
     
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
 MyImage = OP(CurrentExperiement.AOSignal,CurrentExperiement.ScanParam,CurrentExperiement.MySimulationBox.z,param.fs_aq,param.c); 
            
save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission','MyImage','R','zR',...
              'DelayLAWS','ActiveLIST','MyImage','Field_Profile','param');

         case 'OS'
%save('C:\Users\mbocoum\Dropbox\PPT - prez\SLIDES_FRANCOIS\scripts\Simulation_fieldOF.mat','AOSignal','ScanParam') 
  MyImage = OS(CurrentExperiement.AOSignal,CurrentExperiement.ScanParam(:,1),...
             CurrentExperiement.ScanParam(:,2),param.df0x,...
             CurrentExperiement.MySimulationBox.z,...
             param.fs_aq,...
             param.c,[param.X0 , param.X1]); 
     
    save(FileName,'x_phantom','y_phantom','z_phantom','MyTansmission',...
                  'DelayLAWS','ActiveLIST','MyImage','Field_Profile','param')            
             
     end
     %%%%%%%%%%%%%%%%%%%%
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End Program - Free memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmpath('..\radon inversion')
% field_end;