%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
field_init(0);

parameters;
IsSaved = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);

% Example : set active elements

 for n_scan = 1:CurrentExperiement.Nscan
if CurrentExperiement.ScanParam(n_scan,1) == 0
Nactive = sum(CurrentExperiement.BoolActiveList(:,n_scan));
Xs        = (0:Nactive-1)*param.width;  

            [~,~,~,EXCITATION] = CalcMatHole(param.f0*1e-6,...
                CurrentExperiement.ScanParam(n_scan,1),...
                CurrentExperiement.ScanParam(n_scan,2),...
                param.nuX0*1e-3,...
                param.nuZ0*1e-3,Xs*1e3,...
                param.fs*1e-6,param.c,...
                param.Bascule ); % Calculer la matrice
end
 end



