classdef AO_FieldBox
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        y
        z
        
        X
        Y
        Z
    end
    
    methods
        function obj = AO_FieldBox(Xrange,Yrange,Zrange,Nx,Ny,Nz)
            
            obj.x = linspace(Xrange(1),Xrange(end),Nx);
            obj.y = linspace(Yrange(1),Yrange(end),Ny);
            obj.z = linspace(Zrange(1),Zrange(end),Nz);
            
          [X,Y,Z] = meshgrid(obj.x,obj.y,obj.z);  
          
          obj.X = X(:);
          obj.Y = Y(:);
          obj.Z = Z(:);

        
        end
        
        function List = Points(obj)
           
            List = [obj.X,obj.Y,obj.Z];

            
            
        end
        
    end
    
end

