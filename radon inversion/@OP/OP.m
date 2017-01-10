classdef OP
    %OP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta
        t
        R  % Input radon transform of the object
        
        w    % dual variable of t in fourier domain
        F_R  % TF (R) with respect to t direction

    end
    
    methods
        function obj = OP(InputImage,theta,t)
            
            if (size(InputImage) == [length(t),length(theta)])
            % checkin that dimension match the input image
            obj.R = InputImage;
            obj.theta = theta;
            obj.t = t;
            
            else
            msg = 'dimension matrix mismatch';
            error(msg) = 'Error occurred.';
            error(msg);
            
            end
            
        end
        
        function [] = Show_R(obj)
           Hf = figure;
           imagesc(obj.theta,obj.t,obj.R)
           xlabel('theta (°)')
           ylabel('t')
           title('Measured Radon Transform')
            
        end
        
        function obj = t_Fourier_R(obj,N)
            % fftshift : swaps the result [0:N/2-1,-N/2:1)] into natural order
            obj.F_R = fftshift(fft(obj.R));
            
            % defining the frequency axis
            f=linspace(-SamplingRate/(2*1540),SamplingRate/(2*1540),length(yt));
            
        end
        
        function [] = Show_TF_R(obj)
           Hf = figure;
           subplot(121)
           imagesc(abs(obj.F_R))
           xlabel('theta (°)')
           ylabel('\omega')
           title('TF (Measured Radon Transform)')
           subplot(122)
           imagesc(angle(obj.F_R))
           xlabel('theta (°)')
           ylabel('\omega')
           title('TF (Measured Radon Transform)')
        end
    
    end

end

