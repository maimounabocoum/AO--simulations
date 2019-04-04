%%% SPPA %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  main  program  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% maimouna bocoum 04-01-2017
clearvars ;

addpath('..\Field_II')
addpath('..\radon inversion')
field_init(0);

parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Start an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CurrentExperiement = Experiment(param);

% initial excitation field :
param.Noc = 250 ;
    t_excitation = (0:1/param.fs:param.Noc*1.5/param.f0);
    excitation   =  sin(2*pi*param.f0*t_excitation).*hanning(length(t_excitation)).^2';
    
%     excitation_env = hilbert(excitation);
%     excitation_env= abs(excitation_env);
% 
    figure;
    plot(t_excitation*1e6,excitation)
    hold on 
    plot(t_excitation*1e6,excitation,'color','red')
    xlabel('time in \mu s')
    ylabel('a.u')
    title('field excitation')
    
    % calcul de 
    P0 = 10e6 ; % max pressure in Pa
    PII = cumtrapz(t_excitation,(P0*excitation).^2/(1.5*1e6)) ; % W/m^2
    PII = PII*1e-4 ; % W/cm^2
    figure; plot(t_excitation*1e6,1e6*PII)
    hold on 
    plot([min(t_excitation*1e6) max(t_excitation*1e6)],1e6*[max(PII)*0.1 max(PII)*0.1],'--')
    hold on
    plot([min(t_excitation*1e6) max(t_excitation*1e6)],1e6*[max(PII)*0.9 max(PII)*0.9],'--')
    
    % get tmin and tmax
    Imin = find(PII>= 0.1*max(PII), 1 ) ;
    Imax = find(PII>= 0.9*max(PII), 1 ) ;
    tmin = t_excitation(Imin)
    tmax = t_excitation(Imax)
    PD = 1.25*(tmax - tmin)
   Isppa = trapz(t_excitation(Imin:Imax),PII(Imin:Imax))/PD 
   xlabel('time \mu s')
   ylabel('\mu J/cm^2')
   title([sprintf('PD %.2f mus', PD*1e6),sprintf('Isppa = %f', Isppa )])
   