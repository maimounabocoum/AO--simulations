classdef LaserBeam
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        k
        z_x
        z_y
        w0_x
        w0_y
        Zr_x
        Zr_y
        
    end
    
    methods
        
        function obj = LaserBeam(param,type)
            
            switch type
                
                case 'gaussian'
                    
            obj.z_x = param.focus(3);
            obj.z_y = param.Rfocus;
            
            w_xi = param.N_elements*param.width/2;
            w_yi = param.element_height/2;
            
            lambda = param.c/param.f0;
            
            obj.k = 2*pi/lambda ; 
            
            delta_X = w_xi^4 - 4*(param.focus(3))^2*lambda^2/pi^2;
            delta_Y = w_yi^4 - 4*(param.Rfocus)^2*lambda^2/pi^2;
            
            obj.w0_x = (w_xi^2 + sqrt(delta_X))/(2);
            obj.w0_y = (w_yi^2 + sqrt(delta_Y))/(2);
            
            obj.Zr_x = pi*(obj.w0_x)^2/lambda;
            obj.Zr_y = pi*(obj.w0_y)^2/lambda;
            
            end
     

        end
        
        function E = GaussianPulse(obj,X,Y,Z)
           
            Z_x = Z - obj.z_x ;
            Z_y = Z - obj.z_y ;
            
            W_x = obj.w0_x*sqrt(1+Z_x .^2/(obj.Zr_x)^2);
            W_y = obj.w0_y*sqrt(1+Z_y.^2/(obj.Zr_y)^2);
            
            XI_x = atan(Z_x/obj.Zr_x);
            R_x = Z_x  + (obj.Zr_x).^2./Z_x ;
            %*exp(-1i*obj.k*Z_x)
            E = (obj.w0_x./W_x).*exp(-X.^2./W_x.^2)...
                .*exp(-1i*(obj.k).*Z_x./(2*R_x)).*exp(1i*XI_x);
            
            
            
        end
        
       
    
    end

end

