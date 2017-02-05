classdef ActuatorProbe
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties   
        rect;
        center; 
        DelayLaw;
    end
    
    properties (Access = private)
   Nactuators
   ActiveList
   end

        
    methods
        
        function obj = ActuatorProbe(Nactuators,Height,Width,no_sub_x,no_sub_y,kerf,ActiveList,Rfocus)
            obj.Nactuators = Nactuators;
            obj.ActiveList = ActiveList;
     
        rect = zeros(length(ActiveList)*no_sub_x*no_sub_y,19);
        center= zeros(Nactuators,3);
               %% absolute center of the probe:
               Xc = (Width + (Nactuators-1)*(kerf+Width))/2;
               %% get center coordinate of all actuator (even if not in ActiveList)
               for i = 1:Nactuators
               center(i,:) = [Width*(i-1)+Width/2,Height/2,0];
               end
               obj.center = Vector_Translation(center,[-Xc,0,0]);
        
               %% initialiaze individual actuators:
                for i = 1:length(ActiveList)                  
                  % creation of a single element with lower left corner located 
                  % at position (0,0,0)
                  Element = SingleElement(Height,Width,no_sub_x,no_sub_y,i);
                  Element = Element_TranslationX(Element,(ActiveList(i)-1)*(kerf+Width));
                  Element = Element_TranslationY(Element,-Height/2);

                  rect([1:no_sub_x*no_sub_y] + no_sub_x*no_sub_y*(i-1),:) = Element;
                end
                  
               %% recentering all elements to (0,0,0):
                       
               rect = Element_TranslationX(rect,-Xc);
               rect = Element_ElevationZ(rect,Rfocus);

                obj.rect = rect;      
                
                
        end
        
        function [] = ShowProbe(obj)
            
            h = figure;
            X = [obj.rect(:,2);obj.rect(:,5);obj.rect(:,8);obj.rect(:,11)];
            Y = [obj.rect(:,3);obj.rect(:,6);obj.rect(:,9);obj.rect(:,12)];
            Z = [obj.rect(:,4);obj.rect(:,7);obj.rect(:,10);obj.rect(:,13)];
            %Element(i,5:7) = [line*Width,(col - 1)*Height,0];
            % Element(i,8:10) = [line*Width,col*Height,0];
            % Element(i,11:13) = [(line - 1)*Width,col*Height,0];
            scatter3(X*1e3,Y*1e3,Z*1e3)
            xlabel('x(mm)')
            ylabel('y(mm)')
            zlabel('z(mm)')
            
        end
          
        function obj = Set_ActuatorDelayLaw(obj,LawSelection,Param,c)
         
              % obj.center : size = Nactuators x 3 (x,y,z coordinates)
              
           switch LawSelection
               % bj.center : corrdinates of actuator center position
               % Param : z focus

               case 'focus'
                   % input parameter = angle of focalisation
               for i = 1:length(obj.ActiveList)
                  Delay(i) = -(1/c)*norm(Param-obj.center(obj.ActiveList(i),:));  
               end 
                  Delay = Delay + norm(Param)/c;
     
               case 'plane'
                   % input parameters = drift angle
                   for i = 1:length(obj.ActiveList)
                   Delay(i) = -(1/c)*(obj.center(obj.ActiveList(i),1)*tan(Param));
                   end
                   
               case 'user'
                   return;
           end
           
           obj.DelayLaw = Delay;
            
        end
 
        %% screening functions %%
        function [] = ShowDelay(obj)
            figure;
            plot(obj.center(obj.ActiveList,1)*1e3,obj.DelayLaw*1e6,'o')
            xlabel('Actuator position z(mm)')
            ylabel('Delay ( \mu s )')
            
        end
        
        function [] = ShowActuatorCenter(obj)
            figure;
            plot(obj.center(:,1)*1e3,'o')
            xlabel('Actuator Index')
            ylabel('position (mm)')
            
        end
        
        
    end
         
    
    
end  
        


        function CenterPosition = Get_ElementCenter(Element)
        
        CenterPosition  = mean(Element(:,17:19));
  
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
        % Z = R + sqrt(R-Y^2)
             
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
    


