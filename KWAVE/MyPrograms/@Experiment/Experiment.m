classdef Experiment
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        
        MyPhantom ;
        MyProbe ;
        MyLaser ;
        MySimulationBox;
        
        % phantom on simulation bow transmission profile
        DiffuseLightTransmission
      
       ScanParam       % Parameter used for the scan
       Nscan           % number of scans performed in the experiement
       BoolActiveList  % active actuator encoded on boolean table, where each column describe 
                       % one scan, and the line index matches the actuator
       AOSignal
        
    end
    
    properties (Access = private)
        % structure containing all the experiement parameters
       param 

    end
    
    methods
        
        % constructor :
        function obj = Experiment(param)
            % param = structure containing all the simulation parameters
                       % for excitation : sparse probe, use with rectangle
            obj.MyPhantom = Phantom( param.phantom.Positions, param.phantom.Sizes, param.phantom.Types );            
            obj.MySimulationBox = AO_FieldBox(param.Xrange,param.Yrange,param.Zrange,param.Nx,param.Ny,param.Nz);

            % IR laser :
            % if no center is specified, the beam center is the 
            if isfield(param,'center')
            obj.MyLaser = LaserBeam(param.w0,param.center);
            else
            obj.MyLaser = LaserBeam(param.w0,obj.MySimulationBox.GetCenter);
            end
            
            obj.MyProbe = ActuatorProbe(param.N_elements,param.element_height,param.width,...
                          param.no_sub_x,param.no_sub_y,param.kerf,1,param.Rfocus);
                      % for waveform :
                      
            obj.param = param;
            
            obj = ConfigureProbeSequence(obj) ;
            
            
        end
        
        function obj = ConfigureProbeSequence(obj)
            
            switch obj.param.FOC_type
                case 'OF'
                    % get the center position X for each scan line               
                       [Scan, obj.ScanParam] = obj.MyProbe.GetIndex( obj.MySimulationBox.x );
                       Scan(isnan(Scan)) = [] ;  
                       % retreive the total number od scans     
                       obj.Nscan = length(Scan);        
                       % initialization : all actuator non-actives
                       obj.BoolActiveList = false(obj.param.N_elements,obj.Nscan);                       
                       % Dimension in index of active probes lenght
                       ActiveWidth = 0.5*obj.param.focus; % as in the experiement, in m
                       ActiveWidth = ceil(ActiveWidth/obj.param.width) ; % convert to index 
                       n_actives = (1:ActiveWidth)  - floor(mean(1:ActiveWidth)) ;
                       % definition of the Translation Matrix :
                       for n_scan = 1:obj.Nscan
                           % to be completed with the width
                           I = Scan(n_scan) + n_actives ;                
                           % removing out of probe elements
                           I(I < 1) = [] ;
                           I(I > obj.param.N_elements) = [] ;   
                           % define transducers as active:
                       obj.BoolActiveList( I , n_scan) = true ; 
                       end                     

                       
               case 'OP'
                    if ~isfield(obj.param,'angles')
                        obj.param.angles = (-45:2:45)*pi/180 ; % default scan angles
                    end
                    
                    % retreive the total number of scans     
                       obj.Nscan = length(obj.param.angles); 
                       obj.ScanParam = obj.param.angles ;
                       obj.BoolActiveList = true(obj.param.N_elements,obj.Nscan);   
                       
                case 'OS'
                    if ~isfield(obj.param,'angles')
                        obj.param.angles = 0 ; % default scan angles
                    end
                    
                    if ~isfield(obj.param,'angles')
                        obj.param.decimation = 2 ; % default decimate
                    end
                    
                    obj.Nscan = length(obj.param.angles)*length(obj.param.decimation); 
                    % we operate the full decimation scan for every
                    % successive angle to scan.
                    [THETA, DEC] = meshgrid(obj.param.angles,obj.param.decimation);
                    obj.ScanParam = [THETA(:),DEC(:)];
                    % initialization : all actuator non-actives
                    obj.BoolActiveList = false(obj.param.N_elements,obj.Nscan); 
                    
                    for i_decimate = 1:length(obj.param.decimation)
                       
                        % selection of column index to modify for a given
                        % decimation
                       I = 1:length(obj.param.decimation):length(obj.param.angles)*length(obj.param.decimation);
                       I = I + (i_decimate-1) ; % index od column with same decimate
                       
                       % set each of those column using decimation map
                       % described by Idecimate 
                       %Idecimate = 1:obj.param.decimation(i_decimate):obj.param.N_elements ;
                       %obj.BoolActiveList( Idecimate , I ) = true ;
                       
                       %Idecimate = 1:obj.param.N_elements ;
                       Imod = mod(1:obj.param.N_elements,2*obj.param.decimation(i_decimate)) ;
                       %Idecimate( Imod < obj.param.decimation(i_decimate) )   = 1 ;
                       %Idecimate( Imod >= obj.param.decimation(i_decimate) )  = 0 ;
                       
                       obj.BoolActiveList( Imod < obj.param.decimation(i_decimate) , I )  = true ;
                       obj.BoolActiveList( Imod >= obj.param.decimation(i_decimate) , I ) = false ;
                    
                    end
 
                    
                    
                
     
            end
            
            % update time of simulation :
           %obj.MySimulationBox.time = 0:(1/obj.param.fs_aq):max(abs(obj.MySimulationBox.z))/(obj.param.c) ;
           
           % data Result Size initialization 
           
           obj.AOSignal = zeros(obj.param.Nz,obj.Nscan) ;
                      
            
        end
            
        function obj = CalculateUSfield(obj,t_excitation,excitation,n_scan)

           FullElementList = 1:obj.param.N_elements ;
           ActiveList =  FullElementList(obj.BoolActiveList(:,n_scan)) ;
           
           % probe strcuture initialization :
           obj.MyProbe = ActuatorProbe(obj.param.N_elements,obj.param.element_height,obj.param.width,...
                          obj.param.no_sub_x,obj.param.no_sub_y,obj.param.kerf,ActiveList,obj.param.Rfocus);

          if obj.param.Activated_FieldII == 1
            % in field, the excitation fieldII should be real

            % define delay law for the probe :
              switch obj.param.FOC_type
                
                 case 'OF' 
                    obj.MyProbe = obj.MyProbe.Set_ActuatorDelayLaw('focus',[obj.ScanParam(n_scan) 0 obj.param.focus],obj.param.c);
                 case 'OP'
                    obj.MyProbe = obj.MyProbe.Set_ActuatorDelayLaw('plane',obj.ScanParam(n_scan),obj.param.c);
                 case 'OS'
                    obj.MyProbe = obj.MyProbe.Set_ActuatorDelayLaw('plane',obj.ScanParam(n_scan,1),obj.param.c);
              end
                    % Initialize home-made probe  :
                    % focus [0 0 0] will be overwritten by the delay law
                    Probe = xdc_rectangles(obj.MyProbe.rect,[10 10 10], [11 0 0]);   

                    % calculate impulse response in FIELD II
                        t_impulseResponse = (0:1/obj.param.fs:2/obj.param.f0);
                        impulse           = sin(2*pi*obj.param.f0*t_impulseResponse);
                        impulse           = impulse.*hanning(length(impulse))'; 
                        xdc_impulse (Probe, impulse);
                    % set excitation in FIELD II:  
                    xdc_excitation (Probe, excitation);
                    % set delay law in FIELD II: 
                    xdc_focus_times(Probe,-1,obj.MyProbe.DelayLaw);
                    % calculate field on MySimulationBox.Points() with FIELD II: 
                    [h,t] = calc_hp(Probe,obj.MySimulationBox.Points());
   
                    % write field results to the current box
                    tmin = t - max(t_excitation)/2;
                    obj.MySimulationBox = obj.MySimulationBox.Get_SimulationResults(tmin,h,obj.param.fs);
                    xdc_free(Probe);
         else
         %%================= implementation without use of FILEDII=====================
         
          switch obj.param.FOC_type
                case 'OF'
                   tmin = (obj.param.Zrange(1)/(obj.param.c)) ;
                   tmax = (max(abs(obj.MySimulationBox.z))/(obj.param.c) + max(t_excitation)) ;
               case 'OP'
                   % width of simulation BOX
                   if obj.ScanParam(n_scan) <= 0
                   DZ0 = (obj.param.Xrange(2))*sin(obj.ScanParam(n_scan)) ;
                   else
                   DZ0 = (obj.param.Xrange(1))*sin(obj.ScanParam(n_scan)) ;
                   end
                   % taking into account additional propagatoion du to tilt scan

                   % angle           
                   tmin = ( obj.param.Zrange(1) )/(obj.param.c) ;
                   %tmin = ( cos(obj.ScanParam(n_scan))*obj.param.Zrange(1) )/(obj.param.c) ;
                   tmax = (max(abs(obj.MySimulationBox.z))/(obj.param.c) + max(t_excitation));
          end
           
           obj.MySimulationBox.time = tmin:(1/obj.param.fs):tmax ;
                                 
           [X,Y,Z] = meshgrid(obj.MySimulationBox.x,obj.MySimulationBox.y,obj.MySimulationBox.z);

           T = (obj.MySimulationBox.time')*ones(1,length(Z(:)));
   
           switch obj.param.FOC_type
                case 'OF'
            [Field,phi]= obj.GaussianPulse(X-obj.ScanParam(n_scan),Y,Z);   
            % F : field to match dimension issued by Field II
            PHI = repmat(phi(:)',length(obj.MySimulationBox.time) , 1 );
            FIELD =  repmat( Field(:)',length(obj.MySimulationBox.time) , 1 ); 
            
            F = interp1(t_excitation,excitation,T - PHI/(2*pi*obj.param.f0)+ max(t_excitation)/2,'linear',0)  ;  

            obj.MySimulationBox.Field = real(F.*FIELD) ;
               case 'OP'
            ZZ = repmat(Z(:)',length(obj.MySimulationBox.time) , 1 );
            XX = repmat(X(:)',length(obj.MySimulationBox.time) , 1 );
                       % F : field to match dimension issued by Field II
            % setting the rotation inveration to center of the simulation
            % box , ie x_c = 0 , and z_c = mean(obj.param.Zrange)
            XI = (XX*sin(obj.ScanParam(n_scan))+ (ZZ- mean(obj.param.Zrange))*cos(obj.ScanParam(n_scan))) ; 
            % rotated variable by angle theta
            
            %figure;
            %imagesc(squeeze( reshape(XI(1,:),[obj.param.Ny,obj.param.Nx,obj.param.Nz]) ) )
            %+ max(t_excitation)/2
            
            %F = imrotate(A,obj.ScanParam(n_scan)) ;
            
            F = interp1(t_excitation,excitation,(T-mean(obj.param.Zrange)/(obj.param.c)) - XI/(obj.param.c),'linear',0)  ;         

            obj.MySimulationBox.Field = real(F) ;  
                   
           end

            %obj.MySimulationBox.Field = evalField('OF'); % OP and OS to be implmented
            
            end
            
            
            
        end 
        %same function to use with parfor   
        function [E,PHI] = GaussianPulse(obj,X,Y,Z)
            
            z_x = obj.param.focus;
            z_y = obj.param.Rfocus;
           
            Z_x = Z - z_x ;
            Z_y = Z - z_y ;
            
            w_xi = obj.param.N_elements*obj.param.width/2;
            w_yi = obj.param.element_height/2;
            
            lambda = obj.param.c/obj.param.f0;
            
             k = 2*pi/lambda ; 
            
            delta_X = w_xi^4 - 4*(z_x)^2*lambda^2/pi^2;
            delta_Y = w_yi^4 - 4*(obj.param.Rfocus)^2*lambda^2/pi^2;
            
            w0_x = (w_xi^2 + sqrt(delta_X))/(2);
            w0_y = (w_yi^2 + sqrt(delta_Y))/(2);
            
            Zr_x = pi*(w0_x)^2/lambda;
            Zr_y = pi*(w0_y)^2/lambda;
               
            W_x = w0_x*sqrt(1+Z_x .^2/(Zr_x)^2);
            W_y = w0_y*sqrt(1+Z_y.^2/(Zr_y)^2);
            
            XI_x = atan(Z_x/Zr_x);
            R_x = Z_x  + (Zr_x).^2./Z_x ;

            E = (w0_x./W_x).*exp(-X.^2./W_x.^2);
                
            % propagation phae + guy phase + sperical wave phase contributions : 
            PHI = k.*Z - XI_x +k*X.^2./(2*R_x);
            
        
            
        end
        
        function obj = EvalPhantom(obj)
            % Evaluate diffuse beam profile on current simulation box :
            % returns a 1x[Nx,Ny,Nz] vector
            DiffuseLight= obj.MyLaser.Eval(obj.MySimulationBox.x,...
                                    obj.MySimulationBox.y,obj.MySimulationBox.z) ;
            % Evaluate absorber position on current simulation box :
            % returns a 1x[Nx,Ny,Nz] vector                     
            Absorbers = obj.MyPhantom.CalculatePhantom(obj.MySimulationBox.x,...
                                                       obj.MySimulationBox.y,obj.MySimulationBox.z) ;
                                                   
            obj.DiffuseLightTransmission = DiffuseLight.*Absorbers ;
          
        end
        
        function [Transmission,R,zR] = ShowPhantom(obj,varargin)
            
            [Nx,Ny,Nz] = obj.MySimulationBox.SizeBox();
            Transmission = squeeze( reshape(obj.DiffuseLightTransmission',[Ny,Nx,Nz]) );
             % check if dimension agree
             if length(obj.MySimulationBox.x) == size(obj.MySimulationBox.x*1e3,2)
                 Transmission = Transmission';
             end
             
             figure;
             if nargin == 2 
                 theta = 180*varargin{1}/pi ;
                 subplot(121)
                 imagesc(obj.MySimulationBox.x*1e3,obj.MySimulationBox.z*1e3,Transmission)
                 xlabel('x(mm)')
                 ylabel('z(mm)')
                 title('diffused light transmission profile')
                 colorbar
                 % interpolate over refined gris before performing
                 % R-transform :
                 dz = (obj.MySimulationBox.z(2) - obj.MySimulationBox.z(1)); % sampling 
                 x = obj.MySimulationBox.x(1):dz:obj.MySimulationBox.x(end);

                 T = interp1(obj.MySimulationBox.x,Transmission',x,'linear',0);
                 [R,zp] = radon(T',theta + 90);
                 zR = mean(obj.MySimulationBox.z) + zp*dz ; % equivalent z sampling 
                 subplot(122)
                 imagesc(theta,zR*1e3,R)
                 ylim(1e3*[min(obj.MySimulationBox.z) max(obj.MySimulationBox.z)])
                 xlabel('\theta(�)')
                 ylabel('z(mm)')
                 title('Radon transform')
                 colorbar
             else
   
             imagesc(obj.MySimulationBox.x*1e3,obj.MySimulationBox.z*1e3,Transmission)
             xlabel('\theta (�)')
             ylabel('z (mm)')
             title('diffused light transmission profile')
             colorbar
             
            end
                
        end
        
        function obj = GetAcquisitionLine(obj,n)

            [Nx,Ny,Nz] = obj.MySimulationBox.SizeBox();
            % full profile calculated here to optimize loop Scan calulation :                   
            LightTransmission = repmat(obj.DiffuseLightTransmission,length(obj.MySimulationBox.time),1)  ;
            
           % MarkedPhotons = (obj.MySimulationBox.Field).^2;
            MarkedPhotons = (obj.MySimulationBox.Field).^2.*LightTransmission ;
            MarkedPhotons = reshape(MarkedPhotons',[Ny,Nx,Nz,length(obj.MySimulationBox.time)]);

%             figure;
%             for i = 1:10:length(obj.MySimulationBox.time)
%             sumY = sum(sum(MarkedPhotons,1),2) ;
%             plot(obj.MySimulationBox.z*1e3,squeeze(sumY(:,:,:,i)))
%             ylim([0 5e-22])
%             title(['c t = ',num2str(obj.MySimulationBox.time(i)*obj.param.c*1e3),'mm'])
%             %imagesc(obj.MySimulationBox.x*1e3,obj.MySimulationBox.z*1e3,squeeze(MarkedPhotons(:,:,:,i))')
%             xlabel('x (mm)')
%             ylabel('z (mm)')
%             colorbar
%             drawnow
%             end

            line = squeeze( sum(sum(sum(MarkedPhotons,1),2),3) );
            % interpolation on simulation box 
            dz_box = obj.MySimulationBox.z(2) - obj.MySimulationBox.z(1) ;
            dz_field = obj.param.c/obj.param.fs ;
            Nint = ceil(dz_box/dz_field); % smoothing parameter which model PhotoDiode?? averaging effect
            % smooth(line,Nint)
            obj.AOSignal(:,n) = interp1((obj.MySimulationBox.time)*obj.param.c,line,obj.MySimulationBox.z,'square',0);
            
%             figure;
%             plot(obj.MySimulationBox.time*obj.param.c*1e3,line,'marker','o')
%             hold on 
%             plot(obj.MySimulationBox.time*obj.param.c*1e3,smooth(line,Nint),'g','linewidth',3)
%             hold on 
%             plot(obj.MySimulationBox.z*1e3,obj.AOSignal(:,n),'--r','linewidth',3)
%             legend({'raw signal','smoothed','interpolated'}) ;
%             title('tagged photons')
%             xlabel('z = c x t (mm)')
%             ylabel('a.u')
            

          
           
            
            
        end
        
        function [] = ShowAcquisitionLine(obj)
               figure;
            switch obj.param.FOC_type
                case 'OF'

            imagesc(obj.ScanParam*1e3,obj.MySimulationBox.z*1e3,obj.AOSignal)
            xlabel('x (mm)')
                case 'OP'
            imagesc(obj.ScanParam*180/pi,obj.MySimulationBox.z*1e3,obj.AOSignal)
            xlabel('angles (�)')
            end
            ylabel('z = ct (mm) ')
            title('\int_{x,y,z} P(x,y,z,t) dxdydz')
            colorbar

            
        end
       
    end
    
end

