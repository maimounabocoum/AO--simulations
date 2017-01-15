classdef OP
    %OP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta
        t  % image dimension in m
        R  % Input radon transform of the object
        
        w    % dual variable of t in fourier domain
        F_R  % TF (R) with respect to t direction

    end
    
    methods
        function obj = OP(InputImage,theta,t,Param,c)
            
            if (size(InputImage) == [length(t),length(theta)])
            % checkin that dimension match the input image
            obj.R = InputImage;
            obj.theta = theta;
            DXsample = c*1/(Param.SamplingRate*1e6) ; 
            obj.t = (1:length(t))*DXsample;
            
            else
            msg = 'dimension matrix mismatch';
            error(msg) = 'Error occurred.';
            error(msg);
            
            end
            
        end
        
        function [] = Show_R(obj)
           Hf = figure;
           imagesc(obj.theta,obj.t*1e3,obj.R)
           xlabel('theta (°)')
           ylabel('t (mm)')
           title('Measured Radon Transform')
            
        end
        
        function L = Get_Lsample(obj)
           L =  max(obj.t) - min(obj.t) ;
        end
        
        function Fm = Fmax(obj)
            dt = obj.t(2) - obj.t(1) ;
            Fm = 1/dt; % in m-1
        end
        
        function obj = t_Fourier_R(obj,N)
            % fftshift : swaps the result [0:N/2-1,-N/2:1)] into natural order
            % performs the zero padding to fit dimmension N
            obj.F_R = fftshift(fft(obj.R,N));
            
            dt = obj.t(2) - obj.t(1) ;
            t = (-N/2:(N/2-1))*dt; % temporal coordinate
            
            % defining the frequency axis
            f = (-N/2:(N/2-1))*(1/dt);
            obj.w = 2*pi*f;
            %linspace(-SamplingRate/(2*1540),SamplingRate/(2*1540),length(obj.t));
            
        end
        
        function [] = Show_TF_R(obj)
           Hf = figure;
           subplot(121)
           imagesc(obj.theta,obj.w,abs(obj.F_R))
           xlabel('theta (°)')
           ylabel('\omega')
           title('TF (Measured Radon Transform)')
           subplot(122)
           imagesc(obj.theta,obj.w,angle(obj.F_R))
           xlabel('theta (°)')
           ylabel('\omega')
           title('TF (Measured Radon Transform)')
        end
    
    end

end

