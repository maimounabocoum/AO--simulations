function [ value ] = fhwm( x ,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x_interpolation = linspace(min(x),max(x),2^14);
y_interpolation = interp1(x,y,x_interpolation,'linear',0);


x_interpolation(y_interpolation < max(y_interpolation)/2 ) = [];
y_interpolation(y_interpolation < max(y_interpolation)/2 ) = [];


value = max(x_interpolation) - min(x_interpolation) ;

% find FWHM :
% figure;
% plot(x,y,'-o')
% hold on 
% plot(x_interpolation,y_interpolation,'red')



end

