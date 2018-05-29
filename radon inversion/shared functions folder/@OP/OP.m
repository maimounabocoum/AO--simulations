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
           figure('position',[520 378 1000 420]);
           
           N_Fc = max( find( abs(obj.f) <= Fc));
           I = (obj.N/2+2):N_Fc ;
           
           subplot(2,2,1)
           imagesc(obj.f( I ),obj.theta*180/pi,abs(obj.F_R( I , :))');
           ylabel('theta (°)')
           xlabel('frequency (m^{-1})')
           title('Fourier transform')
           subplot(2,2,3)
           surfc(obj.l(I(50:end))*1e3,obj.theta*180/pi,abs(obj.F_R( I(50:end) ,:))');
           shading interp
           view(0,90)
           ylabel('theta (°)')
           xlabel('\lambda (mm)')
           title('Fourier transform')
           subplot(2,2,[2 4])
           % find central point for the unwrapping correspondind to theta =
           % 0 , and f = 0 :
           Iy = length(obj.f)/2+1;
           Ix = find(abs(obj.theta) == min(abs(obj.theta)));
           PHASE = UnwrapPHASE(angle(obj.F_R) ,Ix,Iy);
      
           imagesc( obj.theta,obj.f(abs(obj.f) < Fc), PHASE( abs(obj.f)<Fc,:) ) ;   
           colorbar
           xlabel('theta (°)')
           ylabel('frequency (m-1)')
           title('Measured R-Fourier Transform')

        end
                
        function Fm = Fmax(obj)
            dt = obj.t(2) - obj.t(1) ;
            Fm = 1/dt; % in m-1
        end
        
        function [I,z_out] = DataFiltering(obj,Lobject)
            
N       = 2^12;
Lobject = 1*1e-3;
Fc      = 1/Lobject;  % Lobject is the size of the object to detect. Using simple model (sinc function)
                      % we set it to kc = 100/Lobject 
obj = obj.InitializeFourier(N,10*Fc);
%MyImage.Show_R();    % show Radon transform (ie interpolated raw data)
obj.Fmax()        % maximum frequency sampling = 1/dt
                 
%% Nyist principle states the sampling of the object to reconstruct to be such that w > w_max/2 

%%% Fourier Tranform of Radon input image with respect to t
% creating a fourier parameter set using the measure sample size in m

obj.F_R = obj.fourier(obj.R) ;
% MyImage = MyImage.PhaseCorrection(Fc);
% MyImage.Show_F_R(Fc); % Fc : cut-off frequency used for screening


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filtered inverse fourier transform :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter options : 'ram-lak' (default) , 'cosine', 'hamming' , 'hann'
FilterType = 'ram-lak';%'ram-lak' 

filt = FilterRadon(obj.f, obj.N ,FilterType , Fc);
filt = filt(:);
FILTER = filt*ones(1,length(obj.theta));
% hold on
% plot(MyImage.f,filt)

%p = bsxfun(@times, p, H); % faster than for-loop
%  I = MyImage.ifourier(MyImage.F_R);
 I = obj.ifourier(obj.F_R.*FILTER);
 %I = MyImage.ifourier(MyImage.F_R);
% extract image back to initial size :

figure ; imagesc(obj.theta,obj.f,abs(I))
 [I,z_out] = ReduceDataSize( I,'y',obj.t,obj.L);%MyImage.L


        end
        
        

    
    end

end

