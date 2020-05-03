LL=[250E-6, 250E-6, 500E-6, 500E-6, 250E-6, 500E-6];
ww=[5E-6, 10E-6, 5E-6, 10E-6, 5E-6, 10E-6];
dd=[0.2E-6, 0.2E-6, 0.2E-6, 0.2E-6, 0.1E-6, 0.1E-6];


for iii=1:length(LL)
    
    ti = 0;
    tf = 2.50E-9;
    tspan=[ti tf];
    y0=[0; 0; 0];
    
    V=LL(iii)*ww(iii)*dd(iii);
    [T,Y]= ode45(@(t,y) rate_eq(t,y,V),tspan,y0);
    figure(1);
    hold on
    plot(T,Y(:,1));
    title('taþýyýcý yoðunluðu');
    xlabel('zaman');
    ylabel('taþýyýcý');
    hold off
    figure(2);
    hold on
    plot(T,Y(:,2));
    title('foton yoðunluðu ');
    xlabel('zaman');
    ylabel('foton');
    hold off
    figure(3);
    hold on
    plot(T,Y(:,3));
    title('Çýkýþ gücü');
    xlabel('zaman');
    ylabel('Güç');
    hold off
end
figure(1); 
hold on
legend('L=250 \mum,w=5 \mum,d=0,2 \mum','L=250 \mum,w=10 \mum,d=0,2 \mum','L=500 \mum,w=5 \mum,d=0,2 \mum'...
    ,'L=500 \mum,w=10 \mum,d=0,2 \mum','L=250 \mum,w=5 \mum,d=0,1 \mum','L=500 \mum,w=10 \mum,d=0,1 \mum')
hold off

figure(2); 
hold on
legend('L=250 \mum,w=5 \mum,d=0,2 \mum','L=250 \mum,w=10 \mum,d=0,2 \mum','L=500 \mum,w=5 \mum,d=0,2 \mum'...
    ,'L=500 \mum,w=10 \mum,d=0,2 \mum','L=250 \mum,w=5 \mum,d=0,1 \mum','L=500 \mum,w=10 \mum,d=0,1 \mum')
hold off

figure(3); 
hold on
legend('L=250 \mum,w=5 \mum,d=0,2 \mum','L=250 \mum,w=10 \mum,d=0,2 \mum','L=500 \mum,w=5 \mum,d=0,2 \mum'...
    ,'L=500 \mum,w=10 \mum,d=0,2 \mum','L=250 \mum,w=5 \mum,d=0,1 \mum','L=500 \mum,w=10 \mum,d=0,1 \mum')
hold off


function dy = rate_eq(t,y,V)
dy = zeros(3,1);
% d = active layer thickness  [meters]
%d = 2E-7;
I = 0.3;
% te = electron lifetime[seconds]
te = 2.2E-9;
% tp = photon lifetime
tp = 1.6E-12;
% BuyukGama = optical confinement factor
BuyukGama = 0.3;
% a = linewidth broadening factor[GAIN CONSTANT][cm^2 -- m^2]
a = 2.5E-20;
% Betasp = malzemenin spontane emisyon faktoru
Betasp=10^-3;
% B= radiatif rekombinasyon sabiti[CM^3/S -- m^3/s]
B=10E-16;
% c =velocity of light[m/s]
c=3*10^8;
% e = electron charge
e = 1.602E-19;
%ug= grup reflecive index
ug=4;
% n0 = carier density at transparancy[cm^-3 -- m^-3 ]
n0 = 1E24;
% V = volume[m^3]
%V=d*5E-6*250E-6;
% h = Planck Constant
h = 6.626074E-34;
% vg=grup velecity[vg= c/ug]
vg= c/ug;
% am= resonator loss [cm^-1 -- m^-1]
am=45E2;
% nuu = frequancy
nuu =(c/1300E-9);

% Ey(1): carier population
% y(2): photon population
% rate equation for carrier density
% Carriers
dy(1) = (I/e) - y(1)/te - (BuyukGama*c*a*(y(1)/V-n0)/ug)*y(2);
% Photons
dy(2) = (BuyukGama*c*a*(y(1)/V-n0)/ug)*y(2)- y(2)/tp + Betasp*B*y(1)^2/V ;
dy(3) = h*nuu*vg*am*dy(2)/2;%
end