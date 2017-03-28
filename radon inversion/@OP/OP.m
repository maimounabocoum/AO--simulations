classdef OP < TF_t
    %OP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        z
        % image parameters
        theta
        R      % Input radon transform of the object        
        F_R    % Fourier Transform of R with respect to t direction       
        L      % [min max] : dimension of input image in z direction
    end
    
    properties (Access = private)
        SamplingRate         % Samplinf frequency in Hz
        c                    % sound velocity in m/s
        Fc                   % Frequency Cut-Off for screening
        Linphase_fit      
    end
    
    methods
        function obj = OP(InputImage,theta,t,SamplingRate,c)
            %N: number of points for fourier transform
            obj@TF_t(t);
            
            if (size(InputImage) == [length(t),length(theta)])
            % checkin that dimension match the input image
            obj.SamplingRate = SamplingRate ;
            obj.z = t; % longitudinal index for reconstruction box
            obj.c = c ;
            obj.L = [min(t),max(t)] ;
            obj.R = InputImage;
            obj.theta = theta;            
    
            else
            msg = 'dimension matrix mismatch';
            error(msg) = 'Error occurred.';
            error(msg);
            
            end
            
        end
        
        function obj = InitializeFourier(obj,varargin)

            if nargin == 2
                N = varargin{1};
                t = obj.t ;    
                % udate fourier parameters
                DXsample = obj.c*1/(obj.SamplingRate) ; 
                obj = obj.Initialize(N,1/DXsample);
                % interpolated trace on fourier param
                obj.R = interp1(t,obj.R,obj.t,'linear',0);
            elseif nargin == 3
                N = varargin{1};
                Fc = varargin{2};
                t = obj.t ;  
                obj = obj.Initialize(N,Fc);        
                % interpolated trace on fourier param
                obj.R = interp1(t,obj.R,obj.t,'linear',0);
                
            end
        
        end
        
        function obj = PhaseCorrection(obj,Fc)
            % linear fit of the phase for each angle :
           Iy = obj.N/2+1;
           %get phase for position theta  = 0
           N_theta0 = find(obj.theta == 0);
           phase0  = repmat( unwrap(angle( obj.F_R( : ,N_theta0  ) )), 1 , length(obj.theta) );
           obj.F_R = abs(obj.F_R).*exp(1i*phase0) ;

            
        end
        
        function [] = Show_R(obj)
           Hf = figure;
           imagesc(obj.theta*180/pi,obj.t*1e3,obj.R)
           xlabel('theta (°)')
           ylabel('t (mm)')
           ylim(1e3*obj.L)
           title('Measured Radon Transform')
           colorbar
            
        end
        
        function [] = Show_F_R(obj,Fc)
           figure;
           
           subplot(121)
           imagesc(obj.theta*180/pi,2*pi*obj.w(2*pi*abs(obj.w) < Fc),abs(obj.F_R(2*pi*abs(obj.w) < Fc,:)));
           xlabel('theta (°)')
           ylabel('frequency (m^{-1})')
           title('Fourier transform')
           
           subplot(122)
           % find central point for the unwrapping correspondind to theta =
           % 0 , and f = 0 :
           Iy = length(obj.f)/2+1;
           Ix = find(abs(obj.theta) == min(abs(obj.theta)));
           PHASE = UnwrapPHASE(angle(obj.F_R) ,Ix,Iy);
           %N_theta0 = find(obj.theta == 0);
           %PHASE = PHASE - repmat(PHASE(:,N_theta0),1,size(PHASE,2));
           
           imagesc( obj.theta,obj.f(abs(obj.f) < Fc), PHASE( abs(obj.f)<Fc,:) ) ;   
           colorbar
           xlabel('theta (°)')
           ylabel('frequency (m-1)')
           title('Measured R-Fourier Transform')
%            figure;
%            subplot(121)
%            PHASE(Iy,1)
%            plot(obj.theta,PHASE(Iy,:),'linewidth',2);
%            title(['unwrapped phase for f = ',num2str(obj.f(Iy))]);
%            subplot(122)
%            for i=1:length(obj.theta)
%            plot(obj.f(abs(obj.f) < Fc),unwrap( angle(obj.F_R(abs(obj.f)<Fc,i) ) ));
%            hold on
%            end
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

