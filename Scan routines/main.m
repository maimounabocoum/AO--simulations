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
 CurrentExperiement.ShowPhantom();

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
     %CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan) ;   % Calculate current AO signal - Photorefractive
     CurrentExperiement = CurrentExperiement.GetAcquisitionLine(n_scan,'holography') ;  % detection choice : photorafractive / holography
     
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
 
 %CurrentExperiement.ShowFFTreconstruction() ;
 %% reconstruction of JM image:
 
Nfft = 2^10;
G = TF2D( Nfft , Nfft , (Nfft-1)*param.nuX0 , (Nfft-1)*param.nuZ0 );
 x_phantom = CurrentExperiement.MySimulationBox.x ;
 z_phantom = CurrentExperiement.MySimulationBox.z - 0.5*max(CurrentExperiement.MySimulationBox.z) ;
 [Xi,Zi] = meshgrid(x_phantom,z_phantom-0*min(z_phantom));
 [X,Z] = meshgrid(G.x,G.z);
 %[MyTansmission,R,zR] = CurrentExperiement.ShowPhantom(param.angles);
ObjectFFT = zeros(Nfft , Nfft);
I_obj = interp2(Xi,Zi,MyTansmission,X,Z,'linear',0);

 I_extract = find(CurrentExperiement.ScanParam(:,1)==0);
 figure;imagesc(CurrentExperiement.AOSignal(:,I_extract + CurrentExperiement.Nscan ))

I_ccd = (2000:4000) ;

% get single NBx Nbz values phase = 0

%I_phase0 = find(CurrentExperiement.ScanParam(:,3)==0.5);

SpectreIN = G.fourier(I_obj);
Spectre= 0*SpectreIN;

 for n_loop = 1:CurrentExperiement.Nscan
 
     myTrace = CurrentExperiement.AOSignal(:,n_loop + CurrentExperiement.Nscan );
     t = CurrentExperiement.AOSignal(:,n_loop);

     Nbx = CurrentExperiement.ScanParam(n_loop,1);
     Nbz = CurrentExperiement.ScanParam(n_loop,2);
%     PHASE = CurrentExperiement.ScanParam(n_loop,3);
     Cnm(n_loop) = sum(myTrace( I_ccd ).*exp(1i*2*pi*Nbz*(param.nuZ0)*(param.c)*t( I_ccd )) );
     
    DecalZ  =   0.39; % ??
    DecalX  =   0; % ??
 
    s =  exp(2i*pi*(DecalZ*Nbz + DecalX*Nbx));

%if Nbz < 10
    ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = -1i*s*Cnm(n_loop);
    ObjectFFT((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );%-s*1i*Cnm(n_loop);
    Spectre((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = SpectreIN((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx);
    Spectre((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( Spectre((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );
    
%end

    %ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) = ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) + s*exp(1i*2*pi*PHASE)*P_tot(n_loop);
    %ObjectFFT((Nfft/2+1)-Nbz,(Nfft/2+1)-Nbx) = conj( ObjectFFT((Nfft/2+1)+Nbz,(Nfft/2+1)+Nbx) );   

 end
 
 %ObjectFFT = abs(Spectre).*exp(1i*angle(ObjectFFT));
 %ObjectFFT = abs(ObjectFFT).*exp(1i*angle(Spectre));
 
Reconstruct = G.ifourier( ObjectFFT );
I_obj_r = G.ifourier( Spectre );
% % I = ifft2(ifftshift(ObjectFFT));
% Reconstruct = Reconstruct - ones(Nfft,1)*Reconstruct(1,:);
% % I = ifftshift(I,2);
figure(2);
subplot(221)
imagesc(G.fx/(param.nuX0),G.fz/(param.nuZ0),abs(ObjectFFT))
axis([-5 5 -23 23])
colorbar
subplot(222)
imagesc(G.x*1e3,G.z*1e3,real(Reconstruct))
ylim([-8 8])
title('reconstructed AO image')
xlabel('x(mm)')
ylabel('z(mm)')
cb = colorbar;
ylabel(cb,'a.u.')
subplot(224)
imagesc(G.x*1e3,G.z*1e3,real(I_obj_r))
ylim([-8 8])
subplot(223)
imagesc(G.fx/(param.nuX0),G.fz/(param.nuZ0),abs(Spectre))
axis([-5 5 -23 23])

%%
figure
spectre1D = Spectre(:,(Nfft/2+1));
spectre1D_simu = ObjectFFT(:,(Nfft/2+1));
%plot( G.fz/(param.nuZ0) , abs(spectre1D)/max(abs(spectre1D)),'o-')
plot( G.fz/(param.nuZ0) , abs(angle(spectre1D.*conj(spectre1D_simu))) ,'o-')
% hold on
% %plot( G.fz/(param.nuZ0) , abs(spectre1D_simu)/max(abs(spectre1D_simu)),'o-')
% plot( G.fz/(param.nuZ0) , unwrap(angle(spectre1D_simu)),'o-')
 xlim([-20 20])
% legend('phantom','simulation')




 %% run this code portion to visualize the field temporal and/or spatial profile
 
 Hf = gcf;      % open a new figure
 n_scan =  1;  % index of the scan - look inside variable for corresponding parameters
 parameters; % script with simulation parameter (to edit befor running the simulation)
 CurrentExperiement = Experiment(param); 
 CurrentExperiement = CurrentExperiement.EvalPhantom();
 CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan)    ;   % Initializes the Probe
 CurrentExperiement = CurrentExperiement.CalculateUSfield(n_scan)   ;   % Calculate the Field Over input BOX
    % % option for screening : XY, Xt , XZt
    % CurrentExperiement.MySimulationBox.ShowMaxField('Xt', Hf);  
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZt',Hf);   
    % CurrentExperiement.MySimulationBox.ShowMaxField('XZ', Hf);
    %
    % CurrentExperiement.MyAO = CurrentExperiement.MyAO.AOsequenceGenerate(param,CurrentExperiement.ScanParam);
     CurrentExperiement.ShowFieldCorrelation('XZ', Hf , 20e-6, 20e-6 ,n_scan); % ('XZ',Hf, startExposure, Exposure time,n_scan)
    % CurrentExperiement.MySimulationBox.ShowMaxField('YZ', Hf);
    
    
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