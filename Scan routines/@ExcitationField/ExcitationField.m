classdef ExcitationField < TF_t
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        t_excitation % time
        EXCITATION % Field for each actuator (1 actuator = 1 line)   
    end
    
    methods
        
        %f0 : defines the central frequency of the accoustique wave:
        function obj =  ExcitationField(varargin)
            
            % Nactive : number of active elements
            % Nactive = sum(obj.BoolActiveList(:,n_scan));
            
            Nactive = varargin{1};       
            FOC_type = varargin{2};
            
               switch FOC_type
             
                case {'OF','OP','OS'} 
                    
                Noc = varargin{3} ;
                f0 = varargin{4} ;
                fs = varargin{5} ;
                
                % creation of fourier structure
            
                
                    % extract OF parameters : 
                obj.t_excitation = (0:1/fs:Noc*1.5/f0);
                excitation   =  sin(2*pi*f0*obj.t_excitation).*hanning(length(obj.t_excitation)).^2'; 
                obj.EXCITATION = repmat(excitation,Nactive,1) ;
                
                case 'JM'
                    
                width = varargin{3} ;
                f0 = varargin{4} ;
                fs = varargin{5} ;
             Xs        = (0:Nactive-1)*width;             % Echelle de graduation en X
            [~,~,~,EXCITATION] = CalcMatHole(f0*1e-6,ScanParam(n_scan,1),ScanParam(n_scan,2),...
                                               nuX0*1e-3,nuZ0*1e-3,Xs*1e3,...
                                               fs*1e-6,c); % Calculer la matrice
            obj.EXCITATION = EXCITATION';
                end         

        end
        
        function [] = ShowField(obj)
           figure;
           imagesc(obj.t_excitation*1e6,1:size(obj.EXCITATION,2),obj.EXCITATION')
           xlabel('element')
           ylabel('time (\mu s)')
           colorbar
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
              Field_interp = interp1(x,F',x_in,'linear', 0 ) ;
              Field_out(:,1,:) = Field_interp ;
              
%               figure;
%               subplot(1,2,1)
%               imagesc(abs(F_TF))
%               subplot(1,2,2)
%               imagesc(angle(Field_interp))

            
            
            
            
            
        end
        
      
    end
    
end

