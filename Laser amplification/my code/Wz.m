function waist = Wz(z_focus,w0,lambda,z_in)
%WZ Summary of this function goes here
%   Detailed explanation goes here

Zr = pi*w0^2/lambda ;

waist = w0*sqrt( 1 + (z_in-z_focus).^2/Zr^2 );


end

