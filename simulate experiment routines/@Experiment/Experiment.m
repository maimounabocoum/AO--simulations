classdef Experiment
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = public)
        MyPhantom ;
        MyProbe ;
        MyLaser ;
        MyUSbeam ;
        MySimulationBox;
    end
    
    properties (Access = private)
        % structure containing all the experiement parameters
       param 
    end
    
    methods
        function obj = Experiment(param)
            % param = structure containing all the simulation parameters
            obj.MyProbe = ActuatorProbe(param.N_elements,param.element_height,param.width,...
                          param.no_sub_x,param.no_sub_y,param.kerf,param.ActiveList,param.Rfocus);
                      
            obj.MyPhantom = Phantom();
            obj.MyLaser = LaserBeam(param,'gaussian');
            obj.MySimulationBox = AO_FieldBox(param.Xrange,param.Yrange,param.Zrange,param.Nx,param.Ny,param.Nz);
            obj.param = param;
        end
        
        function obj = StartExperiement(obj)
            
            %[h,t] = calc_hp(Probe,CurrentExperiement.MySimulationBox.Points());
            %h = h/max(h(:));
            
            % calculate laser on simulation box
            %obj.MySimulationBox = obj.MySimulationBox.Get_SimulationResults(t,h,obj.param.fs);

            obj.MySimulationBox.time = (0:(obj.param.Nz-1))*(1/obj.param.fs);
            
                % adding delay law
%             Noc = 8; % number of optical cycles
%             t_excitation = (0:1/obj.param.fs:Noc*1.5/obj.param.f0);
%             excitation =  sin(2*pi*obj.param.f0*t_excitation);
%             excitation = excitation.*hanning(length(excitation))';
            % 
            [X,Y,Z] = meshgrid(obj.MySimulationBox.x,obj.MySimulationBox.y,obj.MySimulationBox.z);
            %clf;
            for i = 20:1:length(obj.MySimulationBox.time)
            Field = obj.MyLaser.GaussianPulse(X,Y,Z);
            Field = Field.*exp(1i*2*pi*obj.param.f0*(obj.MySimulationBox.time(i)- Z/(obj.param.c))).*...
                   exp(-(obj.MySimulationBox.time(i) - Z/(obj.param.c)).^2/(8/(obj.param.f0))^2);
            %   imagesc(1e3*obj.MySimulationBox.x,1e3*obj.MySimulationBox.z,real(squeeze(Field))')
            %obj.My

            obj.MySimulationBox.Field(i,:) = real(Field(:))';
            %colorbar
            %title(['t = ',num2str(1e6*obj.MySimulationBox.time(i))]);
           
            %drawnow
            
            end
            

            
            %obj.MySimulationBox.Field = evalField('OF'); % OP and OS to be implmented
            
            
        end
       
    end
    
end

