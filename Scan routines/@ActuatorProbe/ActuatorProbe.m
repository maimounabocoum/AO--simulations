classdef ActuatorProbe
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties   
        rect;      %% all points of sampled actuator defined according to FIELDII requirement (p.38 user_guide)
        % 1. index of element tranducer (always=1 here)
        % 2-4.  First corner coordinates
        % 5-7.  Second corner coordinates
        % 8-10. Third corner coordinates
        % 11 apodisation of the element( here always 0)
        
        rectActive;% all active points defined according to FIELDII requirement (p.38 user_guide)
        center;    % list of actuator center coordinates (x,y,z)
        ActiveList % index of active elements.
        DelayLaw;  % emission time of active elements
    end
    
   properties (Access = private)
   Nactuators
   no_sub_x % number of sub division for single actuator along x
   no_sub_y % number of sub division for single actuator along y
   Height   % single actuator height
   Width    % single actuator width
   kerf     % actuator inter-spacing
   Rfocus   % fixed actuator elevation
   end

        
    methods
        

        function obj = ActuatorProbe(Nactuators,Height,Width,no_sub_x,no_sub_y,kerf,Rfocus)
            
            obj.Nactuators = Nactuators   ;
            ActiveList     = 1:Nactuators ;
        
        % construction of probe coordinate matrix (cf p.38)
        rect = zeros(Nactuators*no_sub_x*no_sub_y,19);
        center= zeros(Nactuators,3);
               %% absolute center of the probe:
               Xc = (Width + (Nactuators-1)*(kerf+Width))/2;
               %% get center coordinate of all actuator (even if not in ActiveList)
               for i = 1:Nactuators
               center(i,:) = [(kerf+Width)*(i-1),0,0];
               end
               obj.center = Vector_Translation(center,[-Xc,0,0]);
        
               %% initialiaze individual active actuators:
                for i = 1:Nactuators                 
                  % creation of a single element with lower left corner located 
                  % at position (0,0,0)
                  Element = SingleElement(Height,Width,no_sub_x,no_sub_y,i);
                  % translation along X to reach actuator position
                  Element = Element_TranslationX(Element,(ActiveList(i)-1)*(kerf+Width));
                  % recentering of the probe to user defined center of
                  % element
                  Element = Element_TranslationY(Element,-Height/2+center(i,2));

                  rect((1:no_sub_x*no_sub_y) + no_sub_x*no_sub_y*(i-1),:) = Element;
                end
                  
               %% recentering element coordinate origin to (0,0,0):
                       
               rect = Element_TranslationX(rect,-Xc);
               rect = Element_ElevationZ(rect,Rfocus);

                obj.rect        = rect;   
                obj.Height      = Height ;
                obj.Width       = Width ;
                obj.no_sub_x    = no_sub_x ;
                obj.no_sub_y    = no_sub_y ;
                obj.kerf        = kerf;
                obj.Rfocus      = Rfocus ;
                
                
        end
        
        function obj = Set_ActiveList(obj,ActiveList)
            % check that value are well into 1 and Nelement :
            
             obj.ActiveList = ActiveList;  
             rect           = zeros(length(obj.ActiveList)*obj.no_sub_x*obj.no_sub_y,19);
             center         = zeros(obj.Nactuators,3);
               %% absolute center of the probe:
               Xc = (obj.Width + (obj.Nactuators-1)*(obj.kerf+obj.Width))/2;
               %% get center coordinate of all actuator (even if not in ActiveList)
               for i = 1:obj.Nactuators
               center(i,:) = [obj.Width*(i-1)+obj.Width/2,obj.Height/2,0];
               end
               center = Vector_Translation(center,[-Xc,0,0]);
        
               %% initialize individual active actuators:
                for i = 1:length(obj.ActiveList)               
                  % creation of a single element with lower left corner located 
                  % at position (0,0,0)
                  Element = SingleElement(obj.Height,obj.Width,obj.no_sub_x,obj.no_sub_y,i);
                  Element = Element_TranslationX(Element,(ActiveList(i)-1)*(obj.kerf+obj.Width));
                  Element = Element_TranslationY(Element,-obj.Height/2);

                  rect([1:obj.no_sub_x*obj.no_sub_y] + obj.no_sub_x*obj.no_sub_y*(i-1),:) = Element;
                end
                  
               %% recentering all elements to (0,0,0):
                       
               rect = Element_TranslationX(rect,-Xc);
               rect = Element_ElevationZ(rect,obj.Rfocus);

                obj.rectActive = rect;      
            
        end
          
        function obj = Set_ActuatorDelayLaw(obj,LawSelection,Param,c)
         
              % obj.center : size = Nactuators x 3 (x,y,z coordinates)

           switch LawSelection
               % bj.center : corrdinates of actuator center position
               % Param : z focus

               case 'focus'
                   % input parameter = angle of focalisation
               Delay = zeros(1,obj.Nactuators) ;
               for i = 1:obj.Nactuators
                   %Delay(i) = -(1/c)*norm(Param-obj.center(obj.ActiveList(i),:)) ;  
                  Delay(i) = -(1/c)*norm(Param-obj.center(i,:)) ;        
               end 

               Delay = Delay - min(Delay(obj.ActiveList));
     
               case 'plane'
                   % input parameters = drift angle
%                    if ( sum(obj.ActiveList) == 0 )
%                        Delay = 0*double(obj.ActiveList) ;
%                    else
%                        for i = 1:length(obj.ActiveList)
%                        Delay(i) = -(1/c)*(obj.center(obj.ActiveList(i),1)*sin(Param)) ;
%                        end
                       for i = 1:obj.Nactuators
                       Delay(i) = (1/c)*(obj.center(i,1)*sin(Param)) ;
                       end

                       Delay = Delay - min(Delay);
 
                   
               case 'JM'
                   % input parameters = (fX,fZ)
                   
                   T0 = 1/(c*Param(2)) ;
                   F0 = 1/T0 ;
                   N = floor(T0*Fe);
                   [X,T] = meshgrid(x,t);
                   Mat = (sin(2*pi*F0*(T - alpha*X)) > 0).*sin(2*pi*f0*T) ;
                   
                   for i = 1:length(obj.ActiveList)
                   Delay(i) = -(1/c)*(obj.center(obj.ActiveList(i),1)*tan(Param)) ;
                   end
                   
                   Delay = Delay - min(Delay(obj.ActiveList));
                   
                   
               case 'user'
                   return;

           end
           
           obj.DelayLaw = Delay;

        end
 
        function [n, xn] = GetIndex(obj,x0,x1)
           
            % indexes of position x
            n = 1:length(obj.center(:,1));
            
            n(obj.center(:,1) < x0 | obj.center(:,1) > x1) = [];
            % removing out of range values
            
            %n(isnan(n)) = []; 
            
            xn = obj.center(n,1) ; % position of returned actuator values
 
        end
        
        %% screening functions %%
        function [] = ShowDelay(obj)
            figure;
            plot(obj.center(:,1)*1e3,obj.DelayLaw*1e6,'o')
            xlabel('Actuator position z(mm)')
            ylabel('Delay ( \mu s )')
            
        end
   
        function [] = ShowActuatorCenter(obj)
            
            figure;
            plot(obj.center(:,1)*1e3,'o')
            xlabel('Actuator Index')
            ylabel('position (mm)')
            
        end
        
        function [] = ShowProbe(obj)
            
            % absolute center of the probe:
             %  Xc = (Width + (Nactuators-1)*(kerf+Width))/2;
            
            h = figure;
            ha = axes ;
             
             % obj.rect(:,1)        : index of piezoactuator
             % obj.rect(:,2:4)      : first corner coordinate
             % obj.rect(:,5:7)      : second corner coordinate
             % obj.rect(:,8:10)     : third corner coordinate
             % obj.rect(:,11:13)    : fourth corner coordinate
             % obj.rect(:,15)       : width of rectangle directtion x = (Ex : 0.2mm/(no_sub_x))
             % obj.rect(:,16)       : height of rectangle directtion y = (Ex : 6mm/(no_sub_y))
             % obj.rect(:,17:19)    : center point of rectangle
             
             
%              X = [obj.rect(:,2),obj.rect(:,5),obj.rect(:,8),obj.rect(:,11),obj.rect(:,17)];
%              Y = [obj.rect(:,3),obj.rect(:,6),obj.rect(:,9),obj.rect(:,12),obj.rect(:,18)] ;
%              Z = [obj.rect(:,4),obj.rect(:,7),obj.rect(:,10),obj.rect(:,13),obj.rect(:,19)];
             X = obj.rect(:,17);
             Y = obj.rect(:,18) ;
             Z = obj.rect(:,19);
             Xc = obj.center(:,1);
             Yc = obj.center(:,2) ;
             Zc = obj.center(:,3);
%              s = ones(1,length(Z(:)));
%              s(obj.ActiveList) = 10;
             scatter3(X(:)*1e3,Y(:)*1e3,Z(:)*1e3,3,Z(:)*1e3);
%              hold on
%              scatter3(X(:)*1e3,repmat(mean(Y(:)),length(X(:)),1)*1e3,Z(:)*1e3,6,Z(:)*1e3);
             hold on
             scatter3(Xc(:)*1e3,Yc(:)*1e3,Zc(:)*1e3,20,Zc(:)*1e3);
             xlabel('X(mm)')
             ylabel('Y(mm)')
             shading interp
             axis equal
             cb = colorbar;
             ylabel(cb,'Z(mm)')
             title('Proce center of rectangles')

            
        end
        
        
    end
         
    
    
end  
        
         
        function Element = SingleElement(Height,Width,no_sub_x,no_sub_y,elemNum)
            
        Height = Height/no_sub_y;
        Width = Width/no_sub_x;
        
            Element = zeros(no_sub_x*no_sub_y,19);
            
            Element(:,1) = elemNum;
            
            for i = 1:no_sub_x*no_sub_y
                
                [col line] = ind2sub([no_sub_y,no_sub_x],i);
                %line: index along x
                %col: index along y
                
            % rectangles coordinates:    
            Element(i,2:4) = [(line - 1)*Width,(col - 1)*Height,0];
            Element(i,5:7) = [line*Width,(col - 1)*Height,0];
            Element(i,8:10) = [line*Width,col*Height,0];
            Element(i,11:13) = [(line - 1)*Width,col*Height,0];
            % apodisation
            Element(i,14) = 1;
            % width
            Element(i,15) = Width;
            %height
            Element(i,16) = Height;
            % center of element
            Element(i,17:19) = [(Element(i,2) + Element(i,5))/2,(Element(i,12) + Element(i,3))/2,0];
            end  
            
        end
        
        function Element = Element_TranslationX(Element,DeltaX)
        
        % addind DeltaX values to all X coordinates:
        
        Element(:,[2,5,8,11,17]) = Element(:,[2,5,8,11,17]) + DeltaX ;
        
        
        end
        
        function Element = Element_TranslationY(Element,DeltaX)
        
        % addind DeltaX values to all X coordinates:
        
        Element(:,[3,6,9,12,18]) = Element(:,[3,6,9,12,18]) + DeltaX ;
        
        
        end
        
        function Vout = Vector_Translation(ListVector,V0)
        % V = [x,y,z] is a vector to which the elements 
        % inside ListVector are translated (ListVector lines = nbre of vectors)
        
        % addind DeltaX values to all X coordinates:

        Vout = ListVector + repmat(V0,[size(ListVector,1),1]);
        
        
        end
        
        function Element = Element_ElevationZ(Element,Rfocus)
        
        % Eq circle in cartesian coordinate : Y^2 + (Z - Rfocus)^2 = R
        % Therfore, if M belongs to circle with coord Y , we have :
        % Z = R +/- sqrt(R-Y^2)
             
        Element(:,[4,7,10,13,19]) = Rfocus - sqrt(Rfocus^2 - Element(:,[3,6,9,12,18]).^2);
        % somehow this is how FIELD II defines the probe, if commented, the
        % results are different ?? 
        Element(:,[4,7,10,13,19]) = Element(:,[4,7,10,13,19]) - max(max(Element(:,[4,7,10,13,19]))) ;
        
        if (2*Rfocus < GetElementHeigth(Element))
            warndlg('Rfocus is smaller than actuator height','!! Warning !!')
        end
        
        end
        
        function H = GetElementHeigth(Element)
        H = max(max(Element(:,[3,6,9,12,18]))) - min(min( Element(:,[3,6,9,12,18]) ));
        end
    


