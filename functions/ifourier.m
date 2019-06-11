function Et=ifourier(Ew,N,T)

    Ew=fftshift(Ew);
    Et=ifft(Ew,N)/T*N;
    Et=ifftshift(Et);
