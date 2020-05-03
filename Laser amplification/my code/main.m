%%% test script %%
clearvars;
%% constants
h = 6.62e-34;% J.s
c = 3e8;
%% parameters sheet
L  = 1e-3;     % crystal length in m
w0  = 100e-6; % active surface 
eta = 0.34; % absorption slope
lambda_p = 532e-9; % pump wavelength
lambda_e = 780e-9; % pump wavelength
nu_p = c/lambda_p;
nu_e = c/lambda_e;
Ep = h*nu_p ;
Ee = h*nu_e ;
sigma_e = 4.1e-19*1e-4; % emission cross section m2 at 790nm
tau = 3.2e-6;    

%% simultion variables
F = TF_t(1024,4e6);
E0s_in    = 0.1e-3; % seed input energy in J
E0p_in    = 10e-3; % pump input energy in J
stau_fwhm = 100e-6; % seed beam
ptau_fwhm = 120e-6; % seed beam
Pulse_in = exp(-log(2)*(2*F.t/stau_fwhm).^6); % seed profile
Pump_in  = exp(-log(2)*(2*F.t/ptau_fwhm).^6); % pup profile
% normalization of input pulse
Pulse_in = E0s_in*Pulse_in/trapz(F.t,Pulse_in); % in W/m^2
Pump_in  = E0p_in*Pump_in/trapz(F.t,Pump_in); % in W
RP = eta*Pump_in/(L*w0^2*Ep) ; % pumping rate

figure(1)
plot(1e6*F.t,Pulse_in)
xlabel('time(\mu s)')
ylabel('seed peak power(W)')
title(['Total energy =',num2str(1e3*trapz(F.t,Pulse_in)),' mJ'])
figure(2)
plot(1e6*F.t,Pump_in)
xlabel('time(\mu s)')
ylabel('pump peak power(W)')
title(['Total =',num2str(trapz(F.t,Pump_in)),' Number/m^3'])
%% difinition of z grid for he simulation
z_grid = linspace(0,L,5000);
dzgrid = z_grid(2) - z_grid(1);

Igrid = repmat(Pulse_in(:)/w0^2,1,length(z_grid));
figure(3)
imagesc(z_grid*1e3,F.t*1e6,Igrid)
xlabel('mm')
ylabel('\mu s')

DeltaN = 0*Igrid;

Is = (Ee)/(tau*sigma_e);

for loop = 2:length(z_grid)
    
     DeltaN(:,loop) = RP(:).*(tau-0.1e-12)./(1 + Igrid(:,loop-1)/Is);
     
     Igrid(:,loop) = Igrid(:,loop-1) + dzgrid*sigma_e*Igrid(:,loop-1).*DeltaN(:,loop) ;
       
%     Igrid(loop,:) = Igrid(loop-1,:) + (F.dt)*c*sigma_e*Igrid(loop-1,:).*DeltaN ...
%                     - c*[0,diff(Igrid(loop-1,:))]; 
    
                        
end

% figure(4)
% imagesc(z_grid*1e3,F.t*1e6,sigma_e*DeltaN)
% title('gain')
% xlabel('z(mm)')
% ylabel('time')
% colorbar

figure(4)
imagesc(z_grid*1e3,F.t*1e6,Igrid./Igrid(:,1))
xlabel('z(mm)')
ylabel('time')
colorbar

figure(1)
hold on
plot(1e6*F.t,Igrid(:,end)*w0^2)
xlabel('time(\mu s)')
ylabel('seed peak power(W)')
title(['Rod length =',num2str(L*1e3),' mm'])
%hold off

%%














