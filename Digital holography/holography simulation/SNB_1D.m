%%% model 1 D
clearvars;


N = 1024;
dx = 5e-6;
F = TF1D(N,1/dx);
Nav = 1:50:300;



for j = 1:length(Nav)

Stemp = zeros(Nav(j),N);
Btemp = zeros(Nav(j),N);
Signal = exp(1i*2*pi*(F.t)/(500e-6));

for i=1:Nav(j)
% Btemp(i,:) = Btemp(i,:) + rand(1,N);
Stemp(i,:) = 50 + real(Signal.*exp(1i*10*pi*rand(1)));
Btemp(i,:) = ( poissrnd(Stemp(i,:)) - real( Stemp(i,:) ) );


%B = mean(Btemp,1);
%S = mean(Stemp+Btemp,1);

%% FFT
S_tf(i,:) = F.fourier( Stemp(i,:) + Btemp(i,:));
B_tf(i,:) = F.fourier( Btemp(i,:) );

end

S_tf = sum(abs(S_tf),1);
B_tf = sum(abs(B_tf),1);

%% filter
fx_c = 1900;
fr = 300 ; % 2400
Filter0     = ((F.f-fx_c).^2 <= (fr)^2);
Imin = find(Filter0, 1 );
Imax = find(Filter0, 1, 'last' );
% Sfiltered        = F.ifourier((S_tf).*Filter0);
% Bfiltered        = F.ifourier((B_tf).*Filter0);

s(j) = trapz(F.f,  S_tf .*Filter0  );
b(j) = trapz(F.f, B_tf.*Filter0 );

fprintf('signal = %s and noise =  %d , S/N = %s \n', s(j), b(j), s(j)/b(j) );

% figure(1)
% plot(F.t,S,'color','red')

figure(2);
hold off
plot(F.f,S_tf,'linewidth',3,'marker','square','color','red')
hold on
plot(F.f,  B_tf ,'linewidth',3,'marker','square','color','blue')
%ylim([0 10e-3])
xlim([500 3000])
title(['spectra for n_{av} = ',num2str(Nav(j))])
line([F.f(Imin) F.f(Imin)],get(gca,'YLim'),'Color',[1 0 0])
line([F.f(Imax) F.f(Imax)],get(gca,'YLim'),'Color',[1 0 0])
drawnow

end

figure(3)
hold on
semilogy(Nav,s./b,'-o')
xlabel('averaging number Nav')
ylabel('SNR')
grid on




