classdef TF2D
    %TF2D Summarz of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        N       % number of point for fourier transform
        xRange  
        x       % position in m
        fx       %frequencz varaiable
        kx       %frequencz variable in rad
        dx
        dfx
        zRange  
        z       % position in m
        fz       %frequencz varaiable
        kz       %frequencz variable in rad
        dz
        dfz
    end
    
    methods
        
        function obj = TF2D(N,Fmax_x,Fmax_z)
            % Fmax corresponds to the mas sampling of image
            
            obj.N = N;
            obj.dx = 1/Fmax_x; % in m
            obj.dz = 1/Fmax_z; % in m
            obj.x = (-N/2:1:N/2-1)*obj.dx;
            obj.z = (-N/2:1:N/2-1)*obj.dz;
            obj.xRange = (N-1)*obj.dx;
            obj.zRange = (N-1)*obj.dz;
            obj.dfx = 1/obj.xRange;
            obj.dfz = 1/obj.zRange;
            obj.fx = (-N/2:1:N/2-1)*obj.dfx;
            obj.fz = (-N/2:1:N/2-1)*obj.dfz;
            obj.kx = 2*pi*obj.fx;
            obj.kz = 2*pi*obj.fz;
            
        end
        
        function Ekxkz = fourier(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxkz=fft2(ifftshift(Exz))*(obj.xRange/obj.N)*(obj.zRange/obj.N) ;
            Ekxkz=fftshift(Ekxkz);
        end
        
         function Exkz = fourierz(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exkz = fft(ifftshift(Exz,1),obj.N,1)*(obj.zRange/obj.N) ;
            Exkz = fftshift(Exkz,1);
        end
        
        function Exz = ifourier(obj, Ekxkz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exz=ifft2(ifftshift(Ekxkz))*(obj.N/obj.xRange)*(obj.N/obj.zRange) ;
            Exz=fftshift(Exz);
        end
        
    end
    
end

