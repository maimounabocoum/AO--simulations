classdef LaserBeam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       w0
    end
    
    methods
        
        function obj = LaserBeam(w0)
            obj.w0 = w0;
        end
        
        function Intensity = Eval(obj,x,y,z)
            
           [X,Y,Z] = meshgrid(x,y,z) ;
           Intensity = exp(-2*(X.^2+Z.^2)/(obj.w0)^2) ;
           
           Intensity = Intensity(:)' ;
           
        end

    end

end

