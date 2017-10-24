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
           
           Intensity = 0*X ;
           
           for i = 1:size(obj.w0,1) 
           
           Intensity = Intensity + ...
               exp(-2*(X-obj.center(i,1)).^2/(obj.w0(i,1))^2 - 2*(Z-obj.center(i,3)).^2/(obj.w0(i,2))^2) ;
           
           
           
           end
           
           Intensity = Intensity(:)' ;
           
        end

    end

end

