%% define colormap
D=[1 1 1;0 0 1;0 1 0;1 1 0;1 0 0;];
F=[0 0.25 0.5 0.75 1];
G=linspace(0,1,256);
cmap=interp1(F,D,G);

D=[0 0 1;0 1 1;1 1 1;1 1 0;1 0 0;];
F=[0 0.25 0.5 0.75 1];
G=linspace(0,1,256);
cmap2=interp1(F,D,G);


%% atomic stream suceptibility:
% Equation (5) , article Poltzik & Khalili 2018, p.3

% Range of analysis
Omega = 10:500 ; % in Hz

%% GW parameters:
kappa_I = 2*pi*500 ;    % LIGO half_bandwith
c       = 3e8 ;
omega0  = 2*pi*c/(1064e-9);
m       = 40 ; %
L       = 4000 ;
Ic      = 840e3;
Xi      = -1./Omega ;        % GW suceptibility

Theta_power = 8*omega0*Ic/(m*c*L);

%% article parameters;
a1 = (kappa_I+1i*Omega)./(kappa_I-1i*Omega) ;


%% define squeezed vacuum state:
x = -5:0.01:5 ;
y = -5:0.01:5 ;
[X,Y] = meshgrid(x,y);
[phi,rho] = cart2pol(X,Y);
r = 0 ; %squeeze factor
A = [X(:)';Y(:)'];

%% evaluation of rotation squuzeing matrix:
theta = pi/4 ;
M = [cos(theta) , -sin(theta) ; sin(theta) , cos(theta)]; % pi/2-rotation

% article Sequino 2021 - to photon formalism
Omega = 10 ;        % frequency in Hz
Omega_sql = 100;    % shift frequency in Hz
gamma_itf  = 1;     % FP linewith
K = (Omega_sql./Omega).^2*gamma_itf^2./(Omega.^2+gamma_itf^2);
M = [1 , 0 ; -K, 1];

B = M*A ;

X = reshape(B(1,:),length(y),length(x));
Y = reshape(B(2,:),length(y),length(x));

% Wigner distibution in phase space (p4/5 article Breintenbach)
S = exp(- X.^2/exp(-2*r) - Y.^2/exp(2*r)) ;

figure;
imagesc(x,y,S)
colormap(cmap)

%% representation of the spin system

Omega = [0:0.1:100,101:1000,1100:100:10000]*(2*pi);

% ideal response using LIGO dimensioning
Omega1 = 3*(2*pi);
Gamma1 = 3*(2*pi);
Response1 = (Gamma1/(2*pi))./( Gamma1^2/4 + (Omega-Omega1).^2) ;

figure(2);
l = area(Omega,Response1,'FaceAlpha',0.1, 'EdgeColor','blue','linewidth',0.01);
set(gca,'Xscale','log')
grid on
  
%% 
OmegaS = (2*pi)*[0:0.1:100,101:1000,1100:100:10000];

GammaS = ((2*pi*100)^3/(2*pi*500))./OmegaS ;

figure; loglog( OmegaS, GammaS/(2*pi) );
hold on 
loglog( OmegaS, 0*GammaS/(2*pi) + 3/(2*pi) );
xlabel('Lamor frequency \Omega_S/(2 \pi)')
ylabel('Effective Linewidth \Gamma_S/(2*\pi)')


%%




