classdef Phantom
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% original image
        Positions
        SizeTumors
        Types
 
     end
    
    methods
        
        function obj = Phantom(Positions,SizeTumors,Types)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% initialization of input variables %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              obj.Positions =   Positions; 
              obj.SizeTumors  = SizeTumors; 
              obj.Types     =    Types; 

   
        end
        
        function I_abs = CalculatePhantom(obj,x,y,z)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %% initialization of input variables %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [X,Y,Z] = meshgrid(x,y,z);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Phatom construction in m %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           I_abs = ones(size(X,1),size(X,2),size(X,3));
           
               for i_tumor = 1:size(obj.Positions,1)
                   x_c = obj.Positions(i_tumor,1) ;
                   y_c = obj.Positions(i_tumor,2) ;
                   z_c = obj.Positions(i_tumor,3) ;
                   sizeTumor = obj.SizeTumors(i_tumor);
                       switch obj.Types{i_tumor}
                           case 'gaussian'
                               I_abs = -exp(-((X-x_c).^2+(Y-y_c).^2+(Z-z_c).^2)/sizeTumor^2) + I_abs ;
                           case 'square'
                               I_abs = -(abs(X-x_c)<= sizeTumor & (abs(Y-y_c) <= sizeTumor) & (abs(Z-z_c) <= sizeTumor)  ) + I_abs ;
                           case 'fringes'
                               I_abs = I_abs.*sin((pi/sizeTumor).*(x_c*X + y_c*Y + z_c*Z)/norm(obj.Positions(i_tumor,:))).^2 ;
                       end
               end
               
               I_abs = I_abs(:)';
           
           
           
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
        
        function [] = ShowTumor(obj,x,z)
            
           I_abs = obj.CalculatePhantom(x,0,z) ;
           I_abs = squeeze(reshape(I_abs,[1,length(x),length(z)]));
           figure;
           imagesc(x*1e3,z*1e3,I_abs)
           xlabel('x (mm)')
           ylabel('z (mm)')
           title('Absorbption profile of Phantom')
           colorbar
        end

        
    end
    
end

