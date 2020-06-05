clearvars ;
%% cavite en onde plane
L = 5e-6; % cavity length
lambda = 800e-9;
c = 3e8;
omega = 2*pi*c/lambda ;
n = 1 + 1i*10 ;


%% definition of coordinate
z1 = 0:(lambda/100):L ;
z2 = L:(lambda/100):2*L ;

k1 = n*omega/c;
k2 = omega/c;

E1 = exp(1i*k1*z1) - exp(-1i*k1*z1);

B = (exp(1i*k1*L) - exp(-1i*k1*L))/(exp(-1i*k2*L) + exp(1i*k2*L));
E2 = B*( exp(1i*k2*z2) + exp(-1i*k2*z2) );


%% plot field real function
figure(1)
plot(z1,real(E1),'r')
hold on
plot(z2,real(E2))
hold off



