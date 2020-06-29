%%% test script %%
clearvars;
parameters;


%% simultion variables
F = TF_t(1024,0.2e6);
E0s_in    = 10.6e-3;  % seed input energy in J
E0p_in    = 50e-3;   % pump input energy in J
stau_fwhm = 3000e-6; % seed beam
ptau_fwhm = 1000e-6; % pump beam
Pulse_in = exp(-log(2)*(2*F.t/stau_fwhm).^6); % seed profile
Pump_in  = exp(-log(2)*(2*F.t/ptau_fwhm).^6); % pup profile
% normalization of input pulse
Pulse_in = E0s_in*Pulse_in/trapz(F.t,Pulse_in); % in W
Pump_in  = E0p_in*Pump_in/trapz(F.t,Pump_in); % in W
Ipump     = Pump_in/(pi*w0_pump^2);
Ipulse    = Pulse_in/(pi*w0_main^2);

% Rp = eta*Pump_in/(pi*w0^2*c*Ep) ;       % pumping rate nbre/m3/s
% RP = N0*sigma_a*Pump_in/(pi*w0^2*Ep) ; % pumping rate s^{-1}

% sigma_a*tau/(pi*w0^2*Ep)

figure(1)
hold off
plot(1e6*F.t,1e-4*Pulse_in/(pi*w0_main^2))
hold on
yline(1e-4*Is,'-.b');
xlabel('time(\mu s)')
ylabel('seed peak intensity(W/cm^{2})')
title(['Total energy = ',num2str(1e3*trapz(F.t,Pulse_in)),' mJ'])

g0 = sigma_e*sigma_a*tau*N0*max(Ipump)/Ep ;


figure(2)
hold off
plot(1e6*F.t,Pump_in)
hold on
plot(1e6*F.t,Pulse_in)
legend('pump(W)','seed(W)')
xlabel('time(\mu s)')
ylabel('pump power(W)')
title(['g0 = ',num2str(1e-2*g0),' cm^{-1}'])
%% definition of z grid for CW simulation

switch Regime
    case 'CW'
z_grid = linspace(0,L,5000);
dzgrid = z_grid(2) - z_grid(1);

IPULSE = repmat(Ipulse(:),1,length(z_grid));
IPUMP  = repmat(Ipump(:),1,length(z_grid));
DeltaN = 0*IPULSE;
DeltaN0 = 0*IPULSE;

for loop = 2:length(z_grid)
    
    DeltaN0(:,loop) = N0*sigma_a*tau*IPUMP(:,loop-1)./( Ep +sigma_a*tau*IPUMP(:,loop-1) ) ;
    
    Is = Ee*(1+Ep +sigma_a*tau*IPUMP(:,loop-1)/Ep)/(tau*sigma_e);
    
    DeltaN(:,loop) = DeltaN0(:,loop)./(1 + IPULSE(:,loop-1)./Is);
         
     IPULSE(:,loop) = IPULSE(:,loop-1) +...
         dzgrid*sigma_e*IPULSE(:,loop-1).*DeltaN(:,loop) ;
     
     IPUMP(:,loop) = IPUMP(:,loop-1) -...
         dzgrid*sigma_a*IPUMP(:,loop-1).*(N0-DeltaN(:,loop)) ;
     
end

figure(4)
hold off
subplot(2,2,(1:2))
imagesc(z_grid*1e3,F.t*1e6,IPUMP*1e-4)
%imagesc(z_grid*1e3,F.t*1e6,DeltaN)
title('I_{pump} evolution')
xlabel('crystal length (mm)')
ylabel('time (\mu s)')
cb = colorbar ;
ylabel(cb,'[W/cm^2]')
subplot(2,2,3)
plot( z_grid*1e3, 100*(w0_main^2/w0_pump^2)*IPULSE(F.N/2+1 , :)/IPUMP(F.N/2+1,1) )
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('I_{out}/I_{pump}[%]')
subplot(2,2,4)
plot( z_grid*1e3, IPULSE(F.N/2+1 , :)/IPULSE(F.N/2+1 , 1) )
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('I_{out}/I_{in}')

figure(1)
hold on
plot(1e6*F.t,1e-4*IPULSE(:,end))
legend('I_{pulse in}(W/cm^2)','I_{sat}','I_{pulse out}(W/cm^2)')

figure(2)
hold on
plot(1e6*F.t,pi*w0_main^2*IPULSE(:,end))
legend('pump(W)','seed(W)','amplified(W)')

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














