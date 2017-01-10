classdef TumorImage
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% original image
        x % axis in pixels
        y % y axis in pixels
        I % image
        Type
     end
    
    methods
        function obj = TumorImage(WindowSize,Position,SizeTumors,Type)
           obj.x = 1:WindowSize(end);
           obj.x = obj.x - mean(obj.x);
           obj.y = 1:WindowSize(1);
           obj.y = obj.y - mean(obj.y);
           obj.Type = Type;
           
           [X,Y] = meshgrid(obj.x,obj.y);
           obj.I = sparse(length(obj.y),length(obj.x));
           
           for i_tumor = 1:size(Position,1)
               sizeTumor = SizeTumors(i_tumor);
               x_c = Position(i_tumor,2);
               y_c = Position(i_tumor,1);
               
                   switch Type
                       case 'gaussian'
                           obj.I = exp(-((X-x_c).^2+(Y-y_c).^2)/sizeTumor^2)+obj.I;
                       case 'square'
                           obj.I = (abs(X-x_c)<= sizeTumor & (abs(Y-y_c) <= sizeTumor) ) +obj.I;
                   end

           end
           
        end
        
        function [] = ShowTumor(obj)
           figure;
           imagesc(obj.x,obj.y,obj.I)
           xlabel('x')
           ylabel('y')
           title(['tumor'])
           colorbar
        end

        
    end
    
end

