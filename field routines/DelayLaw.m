function Delay_out = DelayLaw( WidthActuator , Nelements,c,focus)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if mod(Nelements,2) == 0
    % case where Nelement is even
    X = WidthActuator *(1:Nelements) - WidthActuator*(Nelements+1)/2 ;

    
    
    
else
    % case where Nelement is even
    
    X = WidthActuator *(1:Nelements) - WidthActuator*Nelements/2 ;
    
end

 % delay in s:
 Delay_out = -(0.3/c)*sqrt(X.^2 + focus^2);
% Delay_out = Delay_out + max(abs(Delay_out));


end

