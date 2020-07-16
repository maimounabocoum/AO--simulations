%%% test script %%
clearvars;
parameters;


%% simultion variables
F = TF_t(1024,6e6);
P0s_in    = 0.150;    % seed input power in W
P0p_in    = 3.2 ;   % pump input power in W
stau_fwhm = 100e-6;   % seed beam
ptau_fwhm = 100e-6;   % pump beam
Pulse_in = exp(-log(2)*4*(2*F.t/stau_fwhm).^6); % seed profile
Pump_in  = exp(-log(2)*4*(2*F.t/ptau_fwhm).^6); % pup profile
% normalization of input pulse
Pulse_in  = P0s_in*Pulse_in; % in W
Pump_in   = P0p_in*Pump_in; % in W
%Pump_in   = E0p_in*Pump_in/trapz(F.t,Pump_in); % in W
Ipump     = Pump_in/(pi*w0_pump^2);
Ipulse    = Pulse_in/(pi*w0_main^2);

% Rp = eta*Pump_in/(pi*w0^2*c*Ep) ;       % pumping rate nbre/m3/s
% RP = N0*sigma_a*Pump_in/(pi*w0^2*Ep) ; % pumping rate s^{-1}

% sigma_a*tau/(pi*w0^2*Ep)

figure(1)
hold off
plot(1e6*F.t,1e-4*Pulse_in/(pi*w0_main^2))
hold on
plot(1e6*F.t,1e-4*repmat(Is,1,F.N),'-.b');
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


z_grid = linspace(0,L,5000);
dzgrid = z_grid(2) - z_grid(1);

IPULSE  = repmat(Ipulse(:),1,length(z_grid));
IPUMP   = repmat(Ipump(:),1,length(z_grid));
DeltaN  = 0*IPULSE;
DeltaN0 = 0*IPULSE;
%%
switch Regime
    case 'CW' % CW case
for loop = 2:length(z_grid)
    
    DeltaN0(:,loop) = N0*sigma_a*tau*IPUMP(:,loop-1)./( Ep +sigma_a*tau*IPUMP(:,loop-1) ) ;
    
    Is = Ee*(1+Ep +sigma_a*tau*IPUMP(:,loop-1)/Ep)/(tau*sigma_e);
    
    DeltaN(:,loop) = DeltaN0(:,loop)./(1 + IPULSE(:,loop-1)./Is);
         
    IPULSE(:,loop) = IPULSE(:,loop-1) +...
         dzgrid*sigma_e*IPULSE(:,loop-1).*DeltaN(:,loop) ;
     
    IPUMP(:,loop) = IPUMP(:,loop-1) -...
         dzgrid*sigma_a*IPUMP(:,loop-1).*(N0-DeltaN(:,loop)) ;
     
end

case 'QCW' % QCW case
for loop = 2:length(z_grid)
    
    % loop on temporal profile
    for loop2 = 2:F.N
       
        DeltaN(loop2,loop) = DeltaN(loop2-1,loop-1) + (F.dt)*( N0*sigma_a*tau*IPUMP(loop2,loop-1)/Ep )...
                             - (F.dt)*( sigma_e*IPULSE(loop2,loop-1)/Ee + sigma_a*IPUMP(loop2,loop-1)/Ep  + 1/tau)*DeltaN(loop2-1,loop-1)   ;
    
    end
    
    % DeltaN0(:,loop) = N0*sigma_a*tau*IPUMP(:,loop-1)./(Ep)- tau*[diff(DeltaN(:,loop-1));0]/(F.dt) ;
    
    % DeltaN(:,loop) = DeltaN0(:,loop)./( tau*sigma_e*IPULSE(:,loop-1)/Ee + 1 + sigma_a*tau*IPUMP(:,loop-1)/Ep );
             
    IPULSE(:,loop) = IPULSE(:,loop-1) + dzgrid*(sigma_e*IPULSE(:,loop-1).*DeltaN(:,loop) - [diff(IPULSE(:,loop-1));0]/(c*F.dt)  ) ;
     
    IPUMP(:,loop) = IPUMP(:,loop-1) - dzgrid*(sigma_a*IPUMP(:,loop-1).*(N0-DeltaN(:,loop)) - [diff(IPUMP(:,loop-1));0]/(c*F.dt) ) ;
     
end


end

%%
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
plot( z_grid*1e3, 100*(w0_main^2/w0_pump^2)*(IPULSE(F.N/2+1 , :)-IPULSE(F.N/2+1 , 1))/IPUMP(F.N/2+1,1) )
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('(I_{out}-I_{in})/I_{pump}[%]')
subplot(2,2,4)
plot( z_grid*1e3, IPULSE(F.N/2+1 , :)/IPULSE(F.N/2+1 , 1) )
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('I_{out}/I_{in}')

figure(1)
hold on
plot(1e6*F.t,1e-4*IPULSE(:,end))
legend('I_{pulse in}(W/cm^2)','I_{sat}','I_{pulse out}(W/cm^2)')

figure(3)
plot(z_grid*1e3,IPUMP(F.N/2+1,:)/IPUMP(F.N/2+1,1));
hold on 
plot(z_grid*1e3,exp(-(z_grid-z_grid(1))*(1e2)*10));
xlabel('mm')
ylabel('a.u')
legend('absorption simulation','th')

figure(2)
hold on
plot(1e6*F.t,pi*w0_main^2*IPULSE(:,end))
legend('pump(W)','seed(W)','amplified(W)')



    












