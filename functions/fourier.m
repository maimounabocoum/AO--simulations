function Ew=fourier(Et,N,T)
    

    t_Et=Et';                %transpose pour avoir la TF sur les lignes
    t_Et=fftshift(t_Et);    %real(F) sera toujours positif pour phi=0
    t_Ew=fft(t_Et,N)*T/N;
    t_Ew=ifftshift(t_Ew);
    Ew=t_Ew';
