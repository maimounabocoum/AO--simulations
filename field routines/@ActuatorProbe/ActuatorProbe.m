classdef ActuatorProbe
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    
        rect;
        

    end

        
    methods
        
        function obj = ActuatorProbe(Nactuators,Height,Width,no_sub_x,no_sub_y,kerf,ActiveList)

        rect = zeros(length(ActiveList)*no_sub_x*no_sub_y,19);
       
        
               %% initialiaze individual actuators:
               for i = 1:length(ActiveList)
                   
                  % creation of a single element with lower left corner located 
                  % at position (0,0,0)
                  Element = SingleElement(Height,Width,no_sub_x,no_sub_y,i);
                  Element = TranslationX(Element,(ActiveList(i)-1)*(kerf+Width));
                  Element = TranslationY(Element,-Height/2);
                  
                  rect(i:(i+no_sub_x*no_sub_y-1),:) = Element;
               end
                  
               % recentering all elements to (0,0,0):
                       % aboslute center :
               Xc = (Width + (Nactuators-1)*(kerf+Width))/2;
               rect = TranslationX(rect,-Xc);
               
        obj.rect = rect;           
        end
               
    end
    
end  
        
        function Element = SingleElement(Height,Width,no_sub_x,no_sub_y,elemNum)
            
            Element = zeros(no_sub_x*no_sub_y,19);
            Element(:,1) = elemNum;
            for i = 1:no_sub_x*no_sub_y
            % rectangles coordinates:    
            Element(i,2:4) = [0,0,0];
            Element(i,5:7) = [Width,0,0];
            Element(i,8:10) = [Width,Height,0];
            Element(i,11:13) = [0,Height,0];
            % apodisation
            Element(i,14) = 1;
            % width
            Element(i,15) = Width;
            %height
            Element(i,16) = Height;
            % center of element
            Element(i,17:19) = [Element(i,2)/2,Element(i,12)/2,0];
            end
            
        end
        
        function Element = TranslationX(Element,DeltaX)
        
        % addind DeltaX values to all X coordinates:
        
        Element(:,[2,5,8,11,17]) = Element(:,[2,5,8,11,17]) + DeltaX ;
        
        
        end
        
        function Element = TranslationY(Element,DeltaX)
        
        % addind DeltaX values to all X coordinates:
        
        Element(:,[3,6,9,12,18]) = Element(:,[3,6,9,12,18]) + DeltaX ;
        
        
        end
    


