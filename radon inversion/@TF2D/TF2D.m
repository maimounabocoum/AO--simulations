classdef TF2D
    %TF2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        N       % number of point for fourier transform
        xRange  
        x       % position in m
        fx       %frequency varaiable
        kx       %frequency variable in rad
        dx
        dfx
        yRange  
        y       % position in m
        fy       %frequency varaiable
        ky       %frequency variable in rad
        dy
        dfy
    end
    
    methods
        
        function obj = TF2D(N,Fmax)
            % Fmax corresponds to the mas sampling of image
            
            obj.N = N;
            obj.dx = 1/Fmax; % in m
            obj.dy = 1/Fmax; % in m
            obj.x = (-N/2:1:N/2-1)*obj.dx;
            obj.y = (-N/2:1:N/2-1)*obj.dy;
            obj.xRange = (N-1)*obj.dx;
            obj.yRange = (N-1)*obj.dy;
            obj.dfx = 1/obj.xRange;
            obj.dfy = 1/obj.yRange;
            obj.fx = (-N/2:1:N/2-1)*obj.dfx;
            obj.fy = (-N/2:1:N/2-1)*obj.dfy;
            obj.kx = 2*pi*obj.fx;
            obj.ky = 2*pi*obj.fy;
            
        end
        
        function Ekxky = fourier(obj, Exy)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxky=fft2(ifftshift(Exy))*(obj.xRange/obj.N)*(obj.yRange/obj.N) ;
            Ekxky=fftshift(Ekxky);
        end
        
        function Exy = ifourier(obj, Ekxky)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exy=ifft2(ifftshift(Ekxky))*(obj.N/obj.xRange)*(obj.N/obj.yRange) ;
            Exy=fftshift(Exy);
        end
        
    end
    
end

