%%% test script %%
clearvars;

param = [0.280,0.5,0.6,0.7,0.8,0.9,1,1.3];          % active surface  ;

myScan = zeros(1,length(param));
I_last = zeros(1,length(param));

 for n_scan = 1:length(param)

    parameters;

%% simultion variables
F = TF_t(2048,5e6);
P0s_in    = param(n_scan);      % 0.28;      % seed input power in W
P0p_in    = 20 ;                % pump input power in W
stau_fwhm = 200e-6;             % seed beam
ptau_fwhm = 170e-6;             % pump beam
Pulse_in = exp(-log(2)*4*(2*(F.t - 88e-6)/stau_fwhm).^200); % seed profile78
Pump_in  = exp(-log(2)*4*(2*F.t/ptau_fwhm).^100);           % pump profile
% difnition of interaction temporal window:
Interact_Window = 1 + 0*F.t;
Interact_Window( Pulse_in < max(Pulse_in)/1e5) = 0 ;

% z grid definition
z_grid = linspace(0,L,500);
x_grid = linspace(0,2*w0_main,100)';

w0_z = Wz( Z_focus , w0_main , lambda_e ,z_grid );
dzgrid = z_grid(2) - z_grid(1);

% normalization of input pulse
Pulse_in  = P0s_in*Pulse_in;    % in W
Pump_in   = P0p_in*Pump_in;     % in W

% intensity profile in x

Xpump_profile  = ( 1/(pi*w0_pump^2) )*exp( -2*(x_grid/w0_pump).^100 )         ;
Xpulse_profile = ( 2/(pi*w0_main^2) )*exp( - 2*(x_grid*(1./w0_main)).^2 )        ;
Xpump_profile  = Xpump_profile/trapz(x_grid,2*pi.*x_grid.*Xpump_profile)      ; % renormalization
Xpulse_profile = Xpulse_profile/trapz(x_grid,2*pi.*x_grid.*Xpulse_profile)      ; % renormalization
%Xpulse_profile = Xpulse_profile./repmat(trapz(x_grid,2*pi.*x_grid.*Xpulse_profile),length(x_grid),1)    ; % renormalization





Ipump     = Xpump_profile*Pump_in ;
Ipulse    = Xpulse_profile*Pulse_in  ;
%Ipulse    = ( 2/(pi*w0_main^2) )*exp( -2*(x_grid/w0_main).^2 )*Pulse_in  ;


% 
 P_pulse_check = trapz( x_grid , 2*pi*repmat(x_grid,1,F.N).*Ipulse ) ;
 P_pump_check = trapz( x_grid , 2*pi*repmat(x_grid,1,F.N).*Ipump ) ;
 figure(1); %hold off ; plot(1e6*F.t, P_pulse_check) ; hold on ; plot(1e6*F.t, Pulse_in);
 plot([-x_grid,x_grid]*1e6,[Xpump_profile,Xpump_profile]) ; hold on ; plot([-x_grid,x_grid]*1e6,[Xpulse_profile(:,end),Xpulse_profile(:,1)]);
 xlabel('r(\mu m)')
 % check profil 
 
 figure(1)
 subplot(2,2,(1:2))
 imagesc( 1e6*F.t , x_grid*1e6 , Ipump*1e-4 )
 cb = colorbar ;
 ylabel(cb,'W/cm^2')
 ylabel('r(\mu m)')
 subplot(2,2,(3:4))
 imagesc( 1e6*F.t , x_grid*1e6 , Ipulse*1e-4 )
 cb = colorbar ;
 ylabel(cb,'W/cm^2')
% hold on
% plot(1e6*F.t,1e-4*repmat(Is,1,F.N),'-.b');
 xlabel('time(\mu s)')
 ylabel('r(\mu m)')
% title(['Total energy = ',num2str(1e3*trapz(F.t,Pulse_in)),' mJ'])
% g0 = sigma_e*sigma_a*tau*N0*max(Ipump)/Ep ;


% figure(2)
% hold off
% plot(1e6*F.t,Pump_in)
% hold on
% plot(1e6*F.t,Pulse_in)
% legend('pump(W)','seed(W)')
% xlabel('time(\mu s)')
% ylabel('peak power(W)')
% title(['g0 = ',num2str(1e-2*g0),' cm^{-1}'])

%% definition of z grid for CW simulation

IPULSE  = repmat(Ipulse,1,1,length(z_grid));
IPUMP   = repmat(Ipump,1,1,length(z_grid));
% check energy conservation:
E_out  = squeeze( trapz( x_grid , 2*pi*repmat(x_grid,1,F.N).*IPULSE ) );



DeltaN  = 0*IPULSE;
DeltaN0 = 0*IPULSE;
%%
switch Regime
    case 'CW' % CW case
for loop = 2:length(z_grid)
        
    DeltaN0(:,:,loop) = N0*sigma_a*tau*IPUMP(:,:,loop-1)./( Ep +sigma_a*tau*IPUMP(:,:,loop-1) ) ;
    
    Is = Ee*( 1 + Ep +sigma_a*tau*IPUMP(:,:,loop-1)/Ep )/(tau*sigma_e);
    
    DeltaN(:,:,loop) = DeltaN0(:,:,loop)./(1 + IPULSE(:,:,loop-1)./Is);
         
    IPULSE(:,:,loop) = IPULSE(:,:,loop-1) +...
         dzgrid*sigma_e*IPULSE(:,:,loop-1).*DeltaN(:,:,loop) ;
     
    IPUMP(:,:,loop) = IPUMP(:,:,loop-1) -...
         dzgrid*sigma_a*IPUMP(:,:,loop-1).*(N0-DeltaN(:,:,loop)) ;
     
     % intensity renormalization:

     IPULSE(:,:,loop) = (w0_z(loop-1)/w0_z(loop))^2*IPULSE(:,:,loop); %
%       P_pulse_check = squeeze(trapz( x_grid , 2*pi*repmat(x_grid,1,1024).*IPULSE  )) ;
%       P_pump_check = squeeze(trapz( x_grid , 2*pi*repmat(x_grid,1,1024).*IPUMP)) ;
%     figure(1)
%     hold on
%     plot(1e6*F.t,P_pulse_check) ;  hold on ; plot(1e6*F.t,P_pump_check) ;
%     colorbar
%     drawnow
     
end

case 'QCW' % QCW case
%%    
for loop = 2:length(z_grid)
    % loop on temporal profile
    a_t = -( (sigma_e/Ee)*IPULSE(:,:,loop-1) + (sigma_a/Ep)*IPUMP(:,:,loop-1) + 1/tau ) ;
    b_t =  N0*(sigma_a/Ep)*IPUMP(:,:,loop-1) ;
    A_t = cumtrapz(F.t,a_t);
    % expmA_t = exp(-A_t) ;
    expmA_t = min(exp(-A_t),realmax('double')); % correction of out of double range error
    
    
    DeltaN(:,:,loop) = cumtrapz( F.t , b_t.*expmA_t )./expmA_t ;

    % DeltaN0(:,loop) = N0*sigma_a*tau*IPUMP(:,loop-1)./( Ep +sigma_a*tau*IPUMP(:,loop-1) ) ;
    % Is = Ee*(1+Ep +sigma_a*tau*IPUMP(:,loop-1)/Ep)/(tau*sigma_e);
    % DeltaN2(:,loop) = DeltaN0(:,loop)./(1 + IPULSE(:,loop-1)./Is);
     
    IPULSE(:,:,loop) = IPULSE(:,:,loop-1) + dzgrid*sigma_e*IPULSE(:,:,loop-1).*DeltaN(:,:,loop) ;
    IPUMP(:,:,loop)  = IPUMP(:,:,loop-1) - dzgrid*sigma_a*IPUMP(:,:,loop-1).*(N0-DeltaN(:,:,loop)) ;
             
    % IPULSE(:,loop) = IPULSE(:,loop-1)  + dzgrid*(sigma_e*IPULSE(:,loop-1).*DeltaN(:,loop)   - [diff(IPULSE(:,loop-1));0]/(c*F.dt)  ) ;     
    % IPUMP(:,loop)  = IPUMP(:,loop-1)   - dzgrid*(sigma_a*IPUMP(:,loop-1).*(N0-DeltaN(:,loop)) - [diff(IPUMP(:,loop-1));0]/(c*F.dt) ) ;
        
    % intensity renomatization:
   IPULSE(:,:,loop) = (w0_z(loop-1)/w0_z(loop))^2*IPULSE(:,:,loop); %
       
    % figure(10)
    % plot(IPULSE(:,loop));
    % drawnow
     
end


end

%%
figure(4)
hold off
subplot(2,2,1)
imagesc(z_grid*1e3,F.t*1e6,squeeze(IPULSE(1,:,:))*1e-4)
%imagesc(z_grid*1e3,F.t*1e6,squeeze(IPUMP(1,:,:))*1e-4)
title('I_{pump} evolution')
xlabel('crystal length (mm)')
ylabel('time (\mu s)')
cb = colorbar ;
ylabel(cb,'[W/cm^2]')
subplot(2,2,2)
imagesc(z_grid*1e3,F.t*1e6,squeeze(IPULSE(:,:,end))*1e-4)
%imagesc(z_grid*1e3,F.t*1e6,squeeze(IPUMP(1,:,:))*1e-4)
title('I_{pump} evolution')
xlabel('crystal length (mm)')
ylabel('time (\mu s)')
cb = colorbar ;
ylabel(cb,'[W/cm^2]')
subplot(2,2,3)

p_pout = squeeze( trapz( x_grid , 2*pi*repmat(x_grid,1,F.N).*IPULSE ) );
p_pump = squeeze( trapz( x_grid , 2*pi*repmat(x_grid,1,F.N).*IPUMP  ) ) ;
E_out  = squeeze( trapz(F.t, p_pout ) );
E_pump = squeeze( trapz(F.t, p_pump ) );

plot( z_grid*1e3, 100*( E_out - E_out(1) )./E_pump(1))
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('(E_{out}- E_{in})/E_{pump}[%]')
subplot(2,2,4)
plot( z_grid*1e3, E_out./E_pump(1) )
title('Rod amplification')
xlabel('crystal length (mm)')
ylabel('P_{out}/P_{in}')

% figure(1)
% hold on
% plot(1e6*F.t,1e-4*IPULSE(:,end))
% legend('I_{pulse in}(W/cm^2)','I_{sat}','I_{pulse out}(W/cm^2)')

% figure(3)
% plot(z_grid*1e3,IPUMP(F.N/2+1,:)/IPUMP(F.N/2+1,1));
% hold on 
% plot(z_grid*1e3,exp(-(z_grid-z_grid(1))*(1e2)*10));
% xlabel('mm')
% ylabel('a.u')
% legend('absorption simulation','th')
% 
% figure(2)
% hold on
% plot(1e6*F.t,pi*(min(w0_z(end),w0_pump))^2*IPULSE(:,end))
% % plot(pi*(min(w0_z(end),w0_pump))^2*IPULSE(:,end))
% legend('pump(W)','seed(W)','amplified(W)')

myScan(n_scan) = p_pout(1296,end);


 end


 figure(3)
 hold on
 plot(param,myScan);
 xlabel('Input Seed Power (W)')
 % ylabel('Extracted Power Ratio (%)')
 % legend('Scan')




