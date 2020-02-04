%% parameter sheet

Fe = 50e6;          % sampling frequency
N = 2^14;
f0 = 6e6;           % US frequency in Hz



F = TF_t(N,Fe);

% structuring frequency 
T0 = 20e-6;             % phase jump duration - caracteristic time
n  = 10;                % harmonic frequency
Np = floor(T0*Fe) ;     % time -> point
N = floor( F.T/T0 );    % number of phase jump
N_c = 10;                % number of period for integration time
nu = 5/( F.T );
t_r = round(rand(1,N));
tau_c = N_c*T0;          % camera integration time

%% difinition de la relation de phase phi(t)
switch Type
    case 'phasejump'
phi = pi + 0*F.t;
As  = 1 + 0*F.t;
for i=0:(N-1)
phi( (i*Np+1):((i+1)*Np)) = phi( (i*Np+1):((i+1)*Np))*t_r(i+1) ;
end

% rajout de la porteuse
% phi = phi + 2*pi*f0*F.t ;

    case 'chirp'


PHI2 = (0.9e6)/tau_c; % 0.185

% phi = 2*pi*f0*F.t + df*(F.t).^2;
 
% periodic ramp

% rajout de la porteuse
% phi = 2*pi*f0*Trep*sawtooth(2*pi*t/Trep) ;
% phi = 2*pi*f0*Trep*sawtooth(2*pi*t/Trep) + 2*pi*df*(Trep*sawtooth(2*pi*t/Trep)).^2;

 phi = 2*pi*PHI2*( T0*( sawtooth(2*pi*F.t/T0) + 1 )/2 ).^2 ;
    case 'JM'
 
 nu = n/T0; % harmonic frequency
 %phi = 0 + 0*F.t; %( 1 + sin(2*pi*F.t/T0) )/2 ;     
 phi = pi*sign( cos(pi*nu*F.t) )/2 ;    
 
% rajout de la porteuse
%  phi = 2*pi*f0*T0*sawtooth(2*pi*F.t/T0) ;

end



















