classdef OP < TF_t
    %OP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta
        R  % Input radon transform of the object        
        F_R  % TF (R) with respect to t direction       
        L
    end
    
    properties (Access = private)
        
        Fc
        Linphase_fit
    end
    
    methods
        function obj = OP(N,InputImage,theta,t,Param,c)
            %N: number of points for fourier transform
            DXsample = c*1/(Param.SamplingRate*1e6) ; 
            obj@TF_t(N,1/DXsample);
            
            if (size(InputImage) == [length(t),length(theta)])
            % checkin that dimension match the input image
            t = (0:(length(t)-1))*DXsample;
            obj.L = max(t);
            obj.R = interp1(t,InputImage,obj.t,'linear',0);
            obj.theta = theta;            
    
            else
            msg = 'dimension matrix mismatch';
            error(msg) = 'Error occurred.';
            error(msg);
            
            end
            
        end
        
        function obj = PhaseCorrection(obj,Fc)
            % linear fit of the phase for each angle :
           Iy = length(obj.f)/2+1;
           %Ix = find(abs(obj.theta) == min(abs(obj.theta)));
           %PHASE = UnwrapPHASE(angle(obj.F_R) ,Ix,Iy);
           phase  = unwrap( angle( obj.F_R( Iy , : ) ) );
        for i=1:length(obj.theta)
           %phase  = unwrap( angle( obj.F_R( abs(obj.f)<Fc , i ) ) );
           %P = polyfit(obj.f(abs(obj.f) < Fc),phase',1);
           %obj.Linphase_fit(:,i) = P(1)*obj.f(abs(obj.f)<Fc)+P(2); 
           %Linphase_fit = smooth(PHASE(Iy,i))*(obj.f)*0.1;
           %Linphase_fit = phase(i) + phase(i)*(obj.f)*2*pi;
           Linphase_fit = (-2+0.00005*obj.theta(i).^2)*abs(obj.f);
        % phase correction 
         obj.F_R(:,i) = obj.F_R(:,i).*exp(-1i*Linphase_fit');
           end 
            
        end
        
        function [] = Show_R(obj)
           Hf = figure;
           imagesc(obj.theta,obj.t*1e3,obj.R)
           xlabel('theta (�)')
           ylabel('t (mm)')
           ylim([0 1e3*obj.L])
           title('Measured Radon Transform')
           colorbar
            
        end
        
        function [] = Show_F_R(obj,Fc)
           Hf = figure;
           subplot(121)
           imagesc(obj.theta,2*pi*obj.w(2*pi*abs(obj.w) < Fc),abs(obj.F_R(2*pi*abs(obj.w) < Fc,:)));
           subplot(122)

           % find central point for the unwrapping correspondind to theta =
           % 0 , and f = 0 :
           Iy = length(obj.f)/2+1;
           Ix = find(abs(obj.theta) == min(abs(obj.theta)));
           PHASE = UnwrapPHASE(angle(obj.F_R) ,Ix,Iy);
           imagesc( obj.theta,obj.f(abs(obj.f) < Fc), PHASE( abs(obj.f)<Fc,:) ) ;   
           colorbar
           xlabel('theta (�)')
           ylabel('frequency (m-1)')
           title('Measured R-Fourier Transform')
           figure;
           subplot(121)
           PHASE(Iy,1)
           plot(obj.theta,PHASE(Iy,:),'linewidth',2);
           title(['unwrapped phase for f = ',num2str(obj.f(Iy))]);
           subplot(122)
           for i=1:length(obj.theta)
           plot(obj.f(abs(obj.f) < Fc),unwrap( angle(obj.F_R(abs(obj.f)<Fc,i) ) ));
           hold on
           end
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
        

    
    end

end

