clearvars
addpath('D:\AO--commons\shared functions folder');
addpath('D:\AO--commons\common subfunctions');


%% generate 1D phase profile

Type = 'chirp'; % 'chirp','periodic'
Parameters;

figure(1)
f0 = 1e6 ; % frequency Delta f en Hz
T_window = 100e-6;
Delfta_f = 1e6;
f = f0 + (Delfta_f/T_window)*(F.t);
E_ref = exp(2*pi*1i*f.*F.t);
H = ones(1,length(F.t));
H(abs(F.t) >= T_window/2) = 0 ;


figure(1); 
subplot(121)
plot(F.t*1e6, unwrap(angle(E_ref.*H.*exp(-1i*2*pi*f0*F.t))) )
%plot(F.t*1e6, real(E_ref.*H))
xlim([-70 70])
xlabel('t( \mu s )')
title('Mudulation')

subplot(122)
S = F.fourier(E_ref.*H);
plot(F.f*1e-6 , abs(S) )
xlim([0 2.5])
xlabel('frequency (MHz)')
title('Mudulation')

DF_eval = (max(f)-min(f))*1e-6


%% modulation d'amplitude
nu0 = 1/T0; % fundamental modulation frequency
%Ar  = 0.5*( 1 + 0*cos( 2*pi*n*nu0*F.t) )  ;
Ar  = cos( 0.5*2*pi*n*nu0*F.t)  ;

%% command AO

par = [0,0.25,0.5,0.75];

%for loop = 1:length(par)

Er = Ar.*exp(1i*( phi ) );
fEr = sign(imag(Er));

z     = 0e-3;  % z point of main tagged photon position
c     = 1540;   % sound velocity in water
delay =  10e-6;     % delayed reference

H  = ( F.t > - tau_c/2  & F.t < tau_c/2  ); % integration window

Es =  interp1( F.t , H.*fEr , F.t -delay + z/c ,'linear',0); % tagged photon at position z / delay [s] = z/v_us 
% 
% figure(100);
% plot(imag(Es)); hold on;
% plot(imag(Er),'r');
% figure(2); hold on ; plot(F.t/T0,real(Er))

% figure(2)
% [h,~,~] = plotyy( F.t/T0 , Ar , F.t/T0 , phi );
% title('reference profile of emission')
% xlabel('time (t/T_0)')
% ylabel(h(1),'amplitude (a.u)')
% ylabel(h(2),'phase (Rad)')

% figure(2);plot(F.t/T0,H);hold on ; plot(F.t/T0,Ar)

% autocorrelation
[acor,lag] = xcorr(Er,Es,'coeff'); % ,
[~,I] = max(abs(acor));
lagDiff = lag(I);
timeDiff = lagDiff/Fe ;
tau = lag/Fe;

% manual evaluation

% T_fwhm = FWHM( abs(acor)/max(abs(acor)) , 1e6*tau );

figure(3);
hold on
plot( tau/T0, abs(acor).^2 )
set(gca,'xtick',[ceil(min(tau/T0)):1:floor(max(tau/T0))])
title('normalized xCorr')
xlabel('\tau / T_0')
ylabel('coefficient')
grid on






% end

% figure(5); hold on ; plot(abs( 1*RR(:,1) + 1i*RR(:,2) - RR(:,3)- 1i*RR(:,4) ) )
% legend(['fwhm = ',num2str(T_fwhm),'\mu s'])

%% imported from mathematica:
%Ith = T0*( 2+cos(2*pi*tau/T0) )/2;

% Ith = T0*( 8*N_c*pi + 4*N_c*pi*cos(2*pi*tau/T0) + 8*sin(N_c*pi) + 4*sin(pi*(N_c-2*tau/T0)) +  sin(2*pi*(N_c + 2*tau/T0)) + 4*sin(pi*(N_c+2*tau/T0))   )/(8*pi);
% hold on
% plot( tau/T0, Ith/(F.dt) )
% legend({'xcorr','mathematica'})

%  hold on
%  Ith = exp(2*1i*pi*PHI2*Trep*(tau-Trep)).*tau.*sinc(2*PHI2*tau.*(tau-Trep)) ...
%      + exp(2*1i*pi*PHI2*Trep*tau).*(Trep-tau).*sinc(2*PHI2*tau.*(tau-Trep)) ;
% %  Ith(tau<0) = 0 ;
% %  Ith(tau>Trep) = 0;
%  plot(1e6*tau,abs(Ith)/max(abs(Ith)))
%  %xlim([0 1e6*Trep])
%  hold off



%%
