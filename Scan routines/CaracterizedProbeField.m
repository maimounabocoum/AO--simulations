%% caracterisation of probe field
%% maimouna bocoum 06-06-2018
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion\shared functions folder')
field_init(0);

parameters;

% overwrite Simulation box initialization : 

    param.Xrange = [-30 30]/1000;  % in m
    param.Yrange = 0/1000;        % [-0.1 0.1]/1000;
    param.Zrange = [4 5]/1000;   % in m

    param.Nx = 150;
    param.Ny = 1;
    % in order to match the number of point in Z direction , and 
    % unshures Nz >=1
    param.Nz = max( 1 , ceil ( param.fs_aq * (abs(param.Zrange(2) - param.Zrange(1)))/(param.c) ) ); % do not edit

% define caracterization box
% x = 0;
% y = 0;
% z = 0;

% initialize parameters
CurrentExperiement = Experiment(param);
t_excitation = (0:1/param.fs:param.Noc*1.5/param.f0);
excitation   =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation)).^2';

Hf = figure;
set(Hf,'WindowStyle','docked'); 
D=[1 1 1;0 0 1;0 1 0;1 1 0;1 0 0;];
F=[0 0.25 0.5 0.75 1];
G=linspace(0,1,256);
cmap=interp1(F,D,G);
colormap(hot)

 for n_scan = 9%:CurrentExperiement.Nscan
 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan);
     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     
     Iout = CurrentExperiement.MySimulationBox.ShowMaxField('Xt',Hf) ;
     subplot(211)
    imagesc(CurrentExperiement.MySimulationBox.x*1e3+19.2,...
           CurrentExperiement.MySimulationBox.time*1e6,Iout*1e12);
                xlabel('x (mm)')
                ylabel('t (\mu s)')
     title(['Nbx = ', num2str(CurrentExperiement.ScanParam(n_scan,2))])           
     cb = colorbar ;
     ylabel(cb,'a.u')
     ylim([2.5 3.5])
     xlim([0 25])
     
 end
 
 %% overwrite Simulation box initialization : 

    param.Zrange = [34 35]/1000;   % in m
    param.Nz = max( 1 , ceil ( param.fs_aq * (abs(param.Zrange(2) - param.Zrange(1)))/(param.c) ) ); % do not edit

%% define caracterization box
% x = 0;
% y = 0;
% z = 0;

% initialize parameters
CurrentExperiement = Experiment(param);
t_excitation = (0:1/param.fs:param.Noc*1.5/param.f0);
excitation   =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation)).^2';

 
     CurrentExperiement = CurrentExperiement.InitializeProbe(n_scan);
     CurrentExperiement = CurrentExperiement.CalculateUSfield(t_excitation,excitation,n_scan);
     
     subplot(212)
     Iout = CurrentExperiement.MySimulationBox.ShowMaxField('Xt') ;
     imagesc(CurrentExperiement.MySimulationBox.x*1e3+19.2,CurrentExperiement.MySimulationBox.time*1e6,Iout*1e12);
                xlabel('x (mm)')
                ylabel('t (\mu s)')                
     cb = colorbar ;
     ylabel(cb,'a.u')
     ylim([22 23])
     xlim([0 25])
     

 set(findall(Hf,'-property','FontSize'),'FontSize',15) 

 saveas(gcf,['C:\Users\mbocoum\Dropbox\self-written documents\acoustic-structured-illumination\images\SimuUSz30',num2str(CurrentExperiement.ScanParam(n_scan,2))],'png')

 
 
 
 
 
 
 
 
 
 
 
 
 