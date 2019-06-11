classdef TF2D
    %TF2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Nx       % number of point for fourier transform
        Ny
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
        
        function obj = TF2D(Nx,Ny,Fex,Fey)
            % Fmax corresponds to the mas sampling of image
            
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.dx = 1/Fex; % in m
            obj.dy = 1/Fey; % in m
            obj.x = (-Nx/2:1:Nx/2-1)*obj.dx;
            obj.y = (-Ny/2:1:Ny/2-1)*obj.dy;
            obj.xRange = (Nx-1)*obj.dx;
            obj.yRange = (Ny-1)*obj.dy;
            obj.dfx = 1/obj.xRange;
            obj.dfy = 1/obj.yRange;
            obj.fx = (-Nx/2:1:Nx/2-1)*obj.dfx;
            obj.fy = (-Ny/2:1:Ny/2-1)*obj.dfy;
            obj.kx = 2*pi*obj.fx;
            obj.ky = 2*pi*obj.fy;
            
        end
        
        function Ekxky = fourier(obj, Exy)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxky = fft2( ifftshift(Exy) , obj.Ny , obj.Nx )*(obj.xRange/obj.Nx)*(obj.yRange/obj.Ny) ;
            %Ekxky = fft2( Exy , obj.N , obj.N )*(obj.xRange/obj.N)*(obj.yRange/obj.N) ;
            Ekxky = fftshift(Ekxky);
        end
        
        function Exy = ifourier(obj, Ekxky)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exy = ifft2(ifftshift(Ekxky),obj.Ny,obj.Nx)*(obj.Nx/obj.xRange)*(obj.Ny/obj.yRange) ;
            Exy = fftshift(Exy);
        end
        
    end
    
end

