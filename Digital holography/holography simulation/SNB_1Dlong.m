%%% model 1 D
clearvars;

N = 512;
Nrep = 5;

for j=1:length(Nrep)
    
dx = 5e-6;
F = TF1D(N*Nrep(j),1/dx);

t = (1:N)*F.dt;
Signal0 = exp(1i*2*pi*t/(500e-6)).*exp(-( (t-mean(t)).^2/(100*F.dt)^2 ) );
Signal = 0*F.t;



    
for i = 1:Nrep(j)

I = ((i-1)*N+1):i*N;
Signal(I) = 50 + real( Signal0.*exp(1i*10*pi*rand(1)) );    
end

Noise = 0.1*( poissrnd(real(Signal)) - real(Signal) );
%Noise = 10*rand(size(Signal));
Signal = Signal + Noise ;

figure(1); 
plot(F.t,Signal,'color','red')


%% fourier transform
SignalTF = F.fourier(Signal);
NoiseTF = F.fourier(Noise);

figure(2)
hold off
plot(F.f,abs(SignalTF).^2,'marker','square','color','red')
hold on
plot(F.f,abs(NoiseTF).^2,'marker','square','color','blue')
xlim([1500 2500])
drawnow 

%% filter
fx_c = 2000;
fr = 400 ; % 2400
Filter0     = ((F.f-fx_c).^2 <= (fr)^2);
Imin = find(Filter0, 1 );
Imax = find(Filter0, 1, 'last' );
% Sfiltered        = F.ifourier((S_tf).*Filter0);
% Bfiltered        = F.ifourier((B_tf).*Filter0);

%s(j) = trapz(F.f,abs( (SignalTF).*Filter0 ).^2 );
s(j) = sum( abs( SignalTF(:).*Filter0(:) ).^2 );
%b(j) = trapz(F.f,abs((NoiseTF).*Filter0).^2);
b(j) = sum( abs( NoiseTF(:).*Filter0(:) ).^2 );
fprintf('signal = %s and noise =  %d , S/N = %s \n',s(j),b(j),s(j)/b(j));


end

figure(3)
hold on
semilogy(Nrep,s./b,'-o')
xlabel('averaging number Nav')
ylabel('SNR')
grid on

