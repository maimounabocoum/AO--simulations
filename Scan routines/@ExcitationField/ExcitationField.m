classdef ExcitationField
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        Probe
    end
    
    properties
        t % time
        omega0
        Excitation % Field for each actuator (1 actuator = 1 line)   
        F0
    end
    
    methods
        
        %f0 : defines the central frequency of the accoustique wave:
        function obj =  ExcitationField(Probe,f0,fs,Noc)
            % input probe intit
            obj.Probe       = Probe ;
            obj.omega0      = 2*pi*f0 ;
            f0              = (obj.omega0)/(2*pi) ;
            obj.omega0      =  2*pi*f0 ;
            obj.t           = (0:1/fs:Noc*1.5/f0);
            obj.Excitation  =  real( exp(1i*obj.omega0*obj.t).*hanning(length(obj.t)).^2' );  

        end
        
        function Field_out = Propagate(obj,x_in,y_in,z_in,c)
            
            Field_out = zeros(length(x_in),length(y_in),length(z_in)) ;
            
            N = 2^10 ;
            txRange = 10*abs(obj.Probe.center(1) - obj.Probe.center(end)) ;
            dx      = txRange/(N-1);
            x   = (-N/2:1:N/2-1)*dx;
            dfx     = 1/txRange ;
            fx      = (-N/2:1:N/2-1)*dfx;
            [Kx,Z] = meshgrid(2*pi*fx,z_in);
            
             % interpolation of inital spatial profile:
             P0 = zeros(1,length(obj.Probe.center));
             P0(obj.Probe.ActiveList) = exp(-1i*2*pi*obj.Probe.DelayLaw*obj.omega0)  ;
             obj.F0 = interp1(obj.Probe.center(:,1),P0,x,'linear',0);
             
             % fourier transform in x :
             F_TF = fft( ifftshift(obj.F0) , N )*txRange/N;
             F_TF = fftshift( F_TF );
             
             F_TF = repmat(F_TF,length(z_in),1);
             % propagation coefficent
  
              F_TF = F_TF.*exp( + 1i*sqrt( (obj.omega0/c)^2 - Kx.^2 ).*Z) ;
              
              % inverse fourier transform
              F = ifft( ifftshift(F_TF,2) , N , 2)*N/txRange;
              F = fftshift(F,2);
             
               % interpolate field over input grid : 
              Field_inter = interp1(x,F',x_in,'linear', 0 ) ;
              
              % so for no implementation in y direction
              Field_out= Field_inter ;
               
%               figure;
%               subplot(1,2,1)
%               imagesc(obj.x,z,abs(F_TF))
%               subplot(1,2,2)
%               imagesc(abs(F))

            
            
            
            
            
        end
        
%         function Ew = fourier(N, Et)
%             %fftshift(Et);    %real(F) sera toujours positif pour phi=0
%             Ew=fft(ifftshift(Et),obj.N)*obj.tRange/obj.N;
%             Ew=fftshift(Ew);
%         end
%         
%         function Et = ifourier(N, Ew)
%             Et=ifft(ifftshift(Ew),N)*obj.N/obj.tRange;
%             Et=fftshift(Et);
%         end
%         
        
    end
    
end

