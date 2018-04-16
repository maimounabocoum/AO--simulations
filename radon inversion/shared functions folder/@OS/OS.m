classdef OS < TF2D
    %OP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % image parameters
        theta
        decimation      % n list 
        ct              % time vector in m
        R               % Input signal interpolated on z grid       
        F_R             % Fourier Transform of R with respect to t direction       
        L               % [min max] : dimension of input image in z direction
    end
    
    properties (Access = private)
        SamplingRate         % Samplinf frequency in Hz
        c                    % sound velocity in m/s
        Fc                   % Frequency Cut-Off for screening
        Linphase_fit      
    end
    
    methods
        
        function obj = OS(InputImage,theta,decimation,df0x,ct,SamplingRate,c)
            %N: number of points for fourier transform
            N = 2^10 ;
            Fmax_x = (N-1)*df0x;
            Fmax_z = 1/(ct(2)-ct(1)) ;
            obj@TF2D(N,Fmax_x,Fmax_z);
            
            if (size(InputImage,1) == length(ct) && size(InputImage,2) == length(theta))
            % checkin that dimension match the input image
            obj.SamplingRate = SamplingRate ;
            obj.ct = ct; % longitudinal index for reconstruction box
            obj.c = c ;
            obj.L = [min(ct),max(ct)] ;
            obj.R = InputImage;
            obj.theta = theta;          
            obj.decimation = decimation; 
            
            % interpolation of raw data into fourier grid :
            obj.R = interp1(ct,obj.R,obj.z,'linear',0);
    
            else
            msg = 'dimension matrix mismatch';
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
                obj.R = interp1(t,obj.R,obj.ct,'linear',0);
            elseif nargin == 3
                N = varargin{1};
                Fc = varargin{2};
                t = obj.ct ;  
                obj = obj.Initialize(N,Fc);        
                % interpolated trace on fourier param
                obj.R = interp1(t,obj.R,obj.ct,'linear',0);
                
            end
        
        end
              
        function [] = Show_R(obj)
           
            % find number of different angular value :
           [Angles,~,Iangles] = unique(obj.theta) ;
           
           Hf = figure;
           
           for iangle=1:length(Angles)
           % for given angle, get proper decimation list :
           Iextract = find(Iangles == iangle) ;
           
           imagesc(obj.decimation(Iextract),obj.z*1e3,obj.R(:,Iextract))
           xlabel('decimation')
           ylabel('ct (mm)')
           ylim(1e3*obj.L)
           title(['Traces for',num2str(Angles(iangle))])
           colorbar
           drawnow
           
           end
            
        end
        
        function [] = Show_F_R(obj)
            
           % find number of different angular value :
           [Angles,~,Iangles] = unique(obj.theta) ;
           
           Hf = figure;
           
           for iangle=1:length(Angles)
           % for given angle, get proper decimation list :
           Iextract = find(Iangles == iangle) ;
                      
           % loop on decimation values to reconstruct fourier composant           
           imagesc(obj.decimation(Iextract),obj.fz*1e-3,abs( obj.F_R(:,Iextract) ) )
           xlabel('decimation')
           ylabel('fz (mm-1)')
           title(['Traces for',num2str(Angles(iangle)*180/pi)])
           colorbar
           drawnow
           
           end
     

        end
        
        function Iout = GetFourierX(obj,Iin,decimation,theta)
            
            % decimation : vector with single element matching decimation
            decimation = unique(decimation);
            [Angles,ia,ib] = unique(theta) ;
            Ntheta = length(Angles);
            I0 = obj.N/2 + 1 ;
            
            % initialization of fourier matrix
            Iout = zeros(obj.N,obj.N,Ntheta) ;
            %
            for n_angle = 1:Ntheta
            
            Iin_angle = Iin(:, find(ib==ia(n_angle)) ) ;   
            
            Iout(:,I0 + decimation,n_angle) = Iin_angle ;
            
            % add conjugate except for 0 order
            % complexe conjuaget
            CONJ = conj(Iin_angle);
            
            for i = 2:length(decimation)
            Iout(:,I0 - decimation(i),n_angle) = [0 ; CONJ(1:end-1,i)] ;
            end
            
            Iout(:,I0,n_angle) = Iout(:,I0,n_angle)/2 ;
            
            end
            
            
        end
        
        function Iout = GetFourier(obj,Iin,decimation,theta)
            
            % decimation : vector with single element matching decimation
            decimation = unique(decimation);
            [Angles,ia,ib] = unique(theta) ;
            Ntheta = length(Angles);
            I0 = obj.N/2 + 1 ;
            
            % initialization of fourier matrix
            Iout = zeros(obj.N,obj.N,Ntheta) ;
            %
            for n_angle = 1:Ntheta
            
            Iin_angle = Iin(:, find(ib==ia(n_angle)) ) ;   
            
            Iout(:,I0 + decimation,n_angle) = Iin_angle ;
            
            % add conjugate except for 0 order
            % complexe conjuaget
            CONJ = conj(flipud(Iin_angle));
            
            for i = 2:length(decimation)
            Iout(:,I0 - decimation(i),n_angle) = [0 ; CONJ(1:end-1,i)] ;
            end
            
            Iout(:,I0,n_angle) = Iout(:,I0,n_angle)/2 ;
            
            end
            
            
        end
        
        function DelayLaw_out = ResizeDelayByAngle(obj,DelayLaw, decimation , theta )
            
            [Decim,iad,ibd] = unique(decimation);
            [Angles,ia,ib] = unique(theta) ;
            Ntheta = length(Angles);
            Ndecimation = length(Decim);
            
            DelayLaw_out = zeros(size(DelayLaw,1),Ntheta,Ndecimation);
            
            for n_dec = 1:Ndecimation

            DelayLaw_out(:,ia,n_dec) = DelayLaw(:,find(n_dec==ibd)) ;

            end
            
        end
        
        function Iout = GetAngles(obj,Iin,decimation,theta)
            
            % decimation : vector with single element matching decimation
            [Decim,iad,ibd] = unique(decimation);
            [Angles,ia,ib] = unique(theta) ;
            Ntheta = length(Angles);
            Ndecimation = length(Decim);
            I0 = obj.N/2 + 1 ;
            
            % initialization of fourier matrix
            Iout = zeros(obj.N,Ntheta,Ndecimation) ;
            %
            for n_dec = 1:Ndecimation

            Iout(:,ia,n_dec) = Iin(:,find(n_dec==ibd)) ;

            end
            
            
        end
        
        function [Iout,theta,decim] = AddSinCos(obj,Iin)
            % Iin : signal 
            % size(Iin,1) = number of lines corresponds to dimention of
            % time vector
            %size(Iin,1) = number of col corresponds to number aquisition
            % this functions takes the four acquisition frame and 
            % compresses it into a single frame according to
            % exp(-ifx x) = cos - i sin
            %   ..        = (h1-h2) - i(h3-h4)
            
           % find number of unique couple values:
           % decimate with be sorted in order
           ScanParam = [obj.decimation(:),obj.theta(:)];
           [Angles,ia,ib] = unique(ScanParam,'rows') ;
           % ia : index of singleton representing group
           % 1:length(thetaUniq) correspond to first 0 decimate
           % ib : index list of goups, same group number is equivalent to
           % same angle and same decimate

           theta = Angles(:,2)  ; % same size as ia 
           thetaUniq = unique(theta,'rows') ;
           
           decim = Angles(:,1)  ; % same size as ia 
           
           % divide the second dimension by the number of phases = 4
           % accounting for single zero order
           % (size(Iin,2)-length(thetaUniq)) : number of decimation
           % + length(thetaUniq) : for 0 order
           Iout = zeros(size(Iin,1),length(thetaUniq) +(size(Iin,2)-length(thetaUniq))/4) ;
           
           Iout(: , 1:length(thetaUniq)) = Iin( : , 1:length(thetaUniq) ) ;
           % starting loop after 0 order
%                        fx = (26.0417)*decim(i_decimate);
%                        Neff = 1/(fx*(0.2*1e-3));
%                        fxeff   = 1/(Neff*1e-3);
%                        

           % length(ia) : number of different groups without repetition
           for i = (length(thetaUniq)+1):length(ia)
               
           Isimilardecimate = sort( find(ib == i) ) ;
           
           % cos = Iin(:,Isimilardecimate(1)) - Iin(:,Isimilardecimate(2))
           % sin = Iin(:,Isimilardecimate(3)) - Iin(:,Isimilardecimate(4))
           
           % sin-cos sequence
           %Iout(:,i) = Iin(:,Isimilardecimate(1)) - 1i*Iin(:,Isimilardecimate(2)) ;
                   
           % own sequence
           
           Iout(:,i) = ( Iin(:,Isimilardecimate(1)) - Iin(:,Isimilardecimate(2)) )...
                     - 1i*( Iin(:,Isimilardecimate(3)) - Iin(:,Isimilardecimate(4)) );
                   
                   
          % Iout(:,i) = hilbert(Iin(:,Isimilardecimate(1)) - Iin(:,Isimilardecimate(2)) );    
            Iout(:,i) = Iout(:,i)/2 ;   
           end       
           
        end
                
        function Fm = Fmax(obj)
            dt = obj.ct(2) - obj.ct(1) ;
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
 [I,z_out] = ReduceDataSize( I,'y',obj.t,obj.L);%MyImage.L


        end
        

    
    end

end

