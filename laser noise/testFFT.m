% test fft relations
N = 50 ;       % number of point for FFT
F = TF_t(N,Fs) ; % Fourier time structure

I = rand(1,N) ; 


%% parseval relationship
I_fft = fft(I) ;

F.dt*sum(I.^2)
(1/(F.N*F.Fs))*sum(abs(I_fft).^2)

%% psd extraction
psd = (1/(F.N*F.Fs))*abs(I_fft(1:F.N0)).^2 ;
psd(2:end-1) = 2*psd(2:end-1) ;

sum(psd)