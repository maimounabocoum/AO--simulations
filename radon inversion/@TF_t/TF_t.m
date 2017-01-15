classdef TF_t
    %TFSpace: Creates an object with all the related variables in the time
    %frequency space
    
    properties
        N       % number of point for fourier transform
        tRange  
        t       % position in m
        f       %frequency varaiable
        w       %frequency variable in rad
        l       %wavelength variable
        dt
        df
        dw
    end
    
    methods 
        function obj = TF_t(N, Fmax)
            % Fmax corresponds to the mas sampling of image

            obj.N = N;
            obj.dt = 1/Fmax; % in m
            obj.t = (-N/2:1:N/2-1)*obj.dt;
            obj.tRange = (N-1)*obj.dt;
            obj.df = 1/obj.tRange;
            obj.f = (-N/2:1:N/2-1)*obj.df;
            obj.w = 2*pi*obj.f;
            obj.l(1:N/2)=-1540./obj.f(1:N/2); % wavelength
            obj.l(N/2+1)=1e15;   % non zero value at origin
            obj.l(N/2+2:N)= 1540./obj.f(N/2+2:N); 
        end
        function Ew = fourier(obj, Et)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ew=fft(ifftshift(Et),obj.N)*obj.tRange/obj.N;
            Ew=fftshift(Ew);
        end
        function Et = ifourier(obj, Ew)
            Et=ifft(ifftshift(Ew),obj.N)*obj.N/obj.tRange;
            Et=fftshift(Et);
        end
    end
    
end

