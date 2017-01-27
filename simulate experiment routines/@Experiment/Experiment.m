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
    
    methods
        function obj = Experiment(param)
            % param = structure containing all the simulation parameters
            obj.MyProbe = ActuatorProbe(param.N_elements,param.element_height,param.width,...
                          param.no_sub_x,param.no_sub_y,param.kerf,param.ActiveList,param.Rfocus);
                      
            obj.MyPhantom = Phantom();
            obj.MyLaser = LaserBeam();
            obj.MySimulationBox = AO_FieldBox(param.Xrange,param.Yrange,param.Zrange,param.Nx,param.Ny,param.Nz);
            
        end
       
    end
    
end

