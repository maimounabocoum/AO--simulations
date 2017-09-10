classdef LaserBeam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       w0
       center 
    end
    
    methods
        
        function obj = LaserBeam(w0,center)
            obj.w0 = w0;
            obj.center = center ;
        end
        
        function Intensity = Eval(obj,x,y,z)
            
           [X,~,Z] = meshgrid(x,y,z) ;
           Intensity = exp(-2*((X-obj.center(1)).^2+(Z-obj.center(3)).^2)/(obj.w0)^2) ;
           
           Intensity = Intensity(:)' ;
           
        end

    end

end

