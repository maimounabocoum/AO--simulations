classdef JM
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nx       % number of point for fourier transform
        Nx0
        Nz
        Nz0
        xRange  
        x        % position in m
        fx       % frequencz varaiable
        kx       % frequencz variable in rad
        dx
        dfx
        zRange  
        z        % position in m
        fz       % frequencz varaiable
        kz       % frequencz variable in rad
        dz
        dfz
    end
    
    methods
        function obj = JM(Nx,Nz,Fmax_x,Fmax_z)
            % Fmax corresponds to the mas sampling of image
            
            obj.Nx = Nx;
            obj.Nz = Nz;
            obj.Nx0 = floor(Nx/2)+1;
            obj.Nz0 = floor(Nz/2)+1;
            obj.dx = 1/Fmax_x; % in m
            obj.dz = 1/Fmax_z; % in m
            obj.x = (-floor(Nx/2):(ceil(Nx/2)-1))*obj.dx;
            obj.z = (-floor(Nz/2):(ceil(Nz/2)-1))*obj.dz;
            obj.xRange = Nx*obj.dx;
            obj.zRange = Nz*obj.dz;
            obj.dfx = 1/obj.xRange;
            obj.dfz = 1/obj.zRange;
            obj.fx = (-floor(Nx/2):(ceil(Nx/2)-1))*obj.dfx;
            obj.fz = (-floor(Nz/2):(ceil(Nz/2)-1))*obj.dfz;
            obj.kx = 2*pi*obj.fx;
            obj.kz = 2*pi*obj.fz;
            
        end
        
        function Ekxkz = fourier(obj, Exz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Ekxkz=fft2(ifftshift(Exz))*(obj.xRange/obj.Nx)*(obj.zRange/obj.Nz) ;
            %Ekxkz=fft2(Exz)*(obj.xRange/obj.Nx)*(obj.zRange/obj.Nz) ;
            Ekxkz=fftshift(Ekxkz);
        end
        
        function Exz = ifourier(obj, Ekxkz)
            %fftshift(Et);    %real(F) sera toujours positif pour phi=0
            Exz=ifft2(ifftshift(Ekxkz))*(obj.Nz/obj.zRange)*(obj.Nx/obj.xRange)  ;
            Exz=fftshift(Exz);
        end
    end
end

