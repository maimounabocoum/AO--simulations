%%
Tmin = 10;
Tmax = 30;
t = 0:0.1:100;
omega1 = 1;
omega2 = 2;
phi = 1.2*pi + 0*t;
phi(t>Tmax) = 0 ;
phi(t>30) = 0 ;

eta = (0.1).*(t>Tmin).*(t<Tmax) ;
E1 = sqrt(1-eta).*exp(1i*omega1*t) + sqrt(eta).*exp(1i*(omega1*t + phi));
E2 = sqrt(1-eta).*exp(1i*(omega2*t) ) + sqrt(eta).*exp(1i*(omega2*t + phi));

figure(1);
hold on
plot(abs(E1+E2).^2)
ylim([0 10])