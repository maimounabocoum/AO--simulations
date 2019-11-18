%% efficicaticite tagg
clearvars;

lambda = 780e-9;
Lambda = 150e-6;
k = 2*pi/lambda ;
K = 2*pi/Lambda;
L = 6e-3;
delta_P = 20e3 ; % Pascal
delta_n = 1.32e-10*delta_P ; % Rad/Pascal

x = k*L*delta_n ;

% Onde progressives
J0  = besselj(0,x)  ;
J1  = besselj(1,x)  ;
J1_ = besselj(-1,x) ;

fprintf('Amplitude oscillation sur collecte J1/J0 : %f percent \n',J1/J0 ) ;
fprintf('Conversion in energy on first order J1 : %f percent \n',100*J1.^2 ) ;

% Onde stationnaires

J0  = besselj(0,x)  ;
J1  = besselj(1,x/2)  ;
J1_ = besselj(-1,x/2) ;

fprintf('Amplitude oscillation sur collecte J1/J0 : %f percent \n',(J1.^2)/(J0.*J0) ) ;

% influence of angle phi
phi = 0:0.001:pi/6;
Lambda = 100e-6;
x = k*L*delta_n*sec(phi).*sinc(pi*L*tan(phi)/Lambda);
x = (2*k*delta_n./(1.33*K*sin(phi))).*sin(0.5*K*L*tan(phi));
figure;
plot(phi*180/pi,besselj(1,x))
xlabel('\phi(°)')
ylabel('J_1(\eta)')



