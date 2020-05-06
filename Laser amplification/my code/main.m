%%% test script %%
clearvars;
parameters;


%% simultion variables
F = TF_t(1024,4e6);
E0s_in    = 0.1e-3; % seed input energy in J
E0p_in    = 10e-3; % pump input energy in J
stau_fwhm = 100e-6; % seed beam
ptau_fwhm = 120e-6; % seed beam
Pulse_in = exp(-log(2)*(2*F.t/stau_fwhm).^6); % seed profile
Pump_in  = exp(-log(2)*(2*F.t/ptau_fwhm).^6); % pup profile
% normalization of input pulse
Pulse_in = E0s_in*Pulse_in/trapz(F.t,Pulse_in); % in W
Pump_in  = E0p_in*Pump_in/trapz(F.t,Pump_in); % in W
%RP = eta*Pump_in/(pi*L*w0^2*Ep) ; % pumping rate nbre/m3/s
RP = sigma_a*Pump_in/(pi*w0^2*Ep) ; % pumping rate s^{-1}



% sigma_a*tau/(pi*w0^2*Ep)

figure(1)
hold off
plot(1e6*F.t,1e-4*Pulse_in/(pi*w0^2))
hold on
yline(1e-4*Is,'-.b');
xlabel('time(\mu s)')
ylabel('seed peak intensity(W/cm^{2})')
title(['Total energy = ',num2str(1e3*trapz(F.t,Pulse_in)),' mJ'])

g0 = sigma_e*max(RP)*tau*N0 ;
g0*1e-2

figure(2)
plot(1e6*F.t,1e-4*Pump_in/(pi*w0^2))
xlabel('time(\mu s)')
ylabel('pump intensity(W/cm^2)')
title(['g0 = ',num2str(1e-2*g0),' cm^{-1}'])
%% definition of z grid for CW simulation
switch Regime
    case 'CW'
z_grid = linspace(0,L,5000);
dzgrid = z_grid(2) - z_grid(1);

Igrid = repmat(Pulse_in(:)/w0^2,1,length(z_grid));
figure(3)
imagesc(z_grid*1e3,F.t*1e6,Igrid)
xlabel('mm')
ylabel('\mu s')

DeltaN = 0*Igrid;

for loop = 2:length(z_grid)
    
     DeltaN(:,loop) = N0*RP(:).*(tau-0.1e-12)./(1 + tau*RP(:) + Igrid(:,loop-1)/Is);
     
     Igrid(:,loop) = Igrid(:,loop-1) + dzgrid*sigma_e*Igrid(:,loop-1).*DeltaN(:,loop) ;
                             
end

figure(4)
hold off
subplot(211)
imagesc(z_grid*1e3,F.t*1e6,Igrid./Igrid(:,1))
xlabel('z(mm)')
ylabel('time')
colorbar
subplot(212)
plot(z_grid*1e3,Igrid(F.N/2+1 , :)/Igrid(F.N/2+1 , 1))
% hold on
% plot(z_grid*1e3,exp(g0*z_grid),'red')
% hold on
% yline(Is/Igrid(F.N/2+1 , 1),'-.b');
xlabel('z(mm)')
ylabel('I_{out}/I_{In}')



    case 'pulsed'
z_grid = linspace(0,L,5000);
dzgrid = z_grid(2) - z_grid(1);
Is = (Ee)/(tau*sigma_e);

% pump population inversion:
[Nu,NL] = PumpInversion(F.t,0,0,RP,tau,0.1e-12);

Igrid = repmat(Pulse_in(:)/w0^2,1,length(z_grid));
figure(3)
imagesc(z_grid*1e3,F.t*1e6,Igrid)
xlabel('mm')
ylabel('\mu s')

DeltaN = 0*Igrid;
DeltaN(:,1) = Nu-NL;

for loop = 2:length(z_grid)
    
     DeltaN(:,loop) = trapz(F.t,-(Igrid(:,loop-1)./(tau*Is)).*(Nu-NL));
     
     Igrid(:,loop) = Igrid(:,loop-1) + dzgrid*sigma_e*Igrid(:,loop-1).*DeltaN(:,loop) ;
                             
end

figure(4)
imagesc(z_grid*1e3,F.t*1e6,DeltaN)
xlabel('z(mm)')
ylabel('time')
colorbar
end














