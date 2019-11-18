%%% model 1 D
clearvars;

N = 512;
Nrep = 2;



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

Noise = poissrnd(real(Signal)) - real(Signal);
%Noise = 10*rand(size(Signal));
Signal = Signal + 0.0000001*Noise ;

figure(1); plot(Signal)

%% Gabor transform
for i_offset = 1:Nrep
wavelet = exp( - ( ( F.t - mean(t)*i_offset).^2/(200*F.dt)^2 ) ) ;
G(i_offset,:) = wavelet.*real(Signal);
%G(i_offset,:) = F.fourier(wavelet.*Signal);
end

%% filter
fx_c = 2000;
fr = 400 ; % 2400
Filter0     = ((F.f-fx_c).^2 <= (fr)^2);
Imin = find(Filter0, 1 );
Imax = find(Filter0, 1, 'last' );

figure(2)
imagesc( F.f(Imin:Imax), 1:size(G,1), abs(G(:,Imin:Imax)))
colorbar

%s(j) = trapz(F.f,abs( (SignalTF).*Filter0 ).^2 );
% s(j) = sum( abs( SignalTF(:).*Filter0(:) ).^2 );
% %b(j) = trapz(F.f,abs((NoiseTF).*Filter0).^2);
% b(j) = sum( abs( NoiseTF(:).*Filter0(:) ).^2 );
% fprintf('signal = %s and noise =  %d , S/N = %s \n',s(j),b(j),s(j)/b(j));


end

% figure(3)
% hold on
% semilogy(Nrep,s./b,'-o')
% xlabel('averaging number Nav')
% ylabel('SNR')
% grid on

