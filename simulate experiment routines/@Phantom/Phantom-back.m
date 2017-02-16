classdef Phantom
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% original image
        x % axis in pixels
        y % y axis in pixels
        I % image
 
     end
    
    methods
        function obj = Phantom(varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% initialization of input variables %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if isempty(find(nargin == 1:4, 1))
              WindowSize = [70,40]*1e-3;
              Position = WindowSize/2;
              SizeTumors = WindowSize/50;
              Type =       'gaussian' ;
          end
           
           if nargin == 1
              WindowSize = varargin{1};
              Position = WindowSize/2;
              SizeTumors = WindowSize/10;
              Type =       'gaussian' ;
           end
            %WindowSize,Position,SizeTumors,Type
            
           if nargin ==4
              WindowSize = varargin{1}; 
              Position =   varargin{2}; 
              SizeTumors = varargin{3}; 
              Type =       varargin{4}; 
           end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Phatom construction in m %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           obj.x = linspace(0,WindowSize(1),500);
           obj.y = linspace(0,WindowSize(end),500);    
           [X,Y] = meshgrid(obj.x,obj.y);
           obj.I = sparse(length(obj.y),length(obj.x));
           
               for i_tumor = 1:size(Position,1)
                   sizeTumor = SizeTumors(i_tumor);
                   x_c = Position(i_tumor,1);
                   y_c = Position(i_tumor,end);

                       switch Type
                           case 'gaussian'
                               obj.I = exp(-((X-x_c).^2+(Y-y_c).^2)/sizeTumor^2)+obj.I;
                           case 'square'
                               obj.I = (abs(X-x_c)<= sizeTumor & (abs(Y-y_c) <= sizeTumor) ) +obj.I;
                       end
               end
           
           
           
        end
        
        function obj = ScanTumor(obj,Angles)
  
            [X,Y] = meshgrid(obj.x,obj.y);
            figure;
            for i=1:length(Angles)
                T = cos(Angles(i)*pi/180)*X + sin(Angles(i)*pi/180)*Y ;
                S = sin(Angles(i)*pi/180)*X - cos(Angles(i)*pi/180)*Y ;
                
                
                imagesc(obj.x*1e3,obj.y*1e3,T)
                xlabel('x (mm)')
                ylabel('y (mm)')
                title(['angle',num2str(Angles(i))])
                drawnow
            end
        end
        
        function [] = ShowTumor(obj)
           figure;
           imagesc(obj.x*1e3,obj.y*1e3,obj.I)
           xlabel('x (mm)')
           ylabel('y (mm) : US probe position')
           title(['Absorbption profile of Phantom'])
           colorbar
        end

        
    end
    
end

