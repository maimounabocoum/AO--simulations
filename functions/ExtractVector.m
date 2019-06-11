function [x_out y_out] = ExtractVector(x_in,y_in,Xmin,Xmax)
%EXTRACTVECTOR Summary of this function goes here
%   Detailed explanation goes here

x_out = x_in(x_in >= Xmin & x_in <= Xmax);


if size(y_in,1)== 1 || size(y_in,2) == 1

y_out = y_in(x_in >= Xmin & x_in <= Xmax);

 else
% 
 y_out = y_in(:,x_in >= Xmin & x_in <= Xmax);
%       
 end


end

