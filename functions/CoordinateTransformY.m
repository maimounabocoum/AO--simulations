function [Xq Yq Fout] = CoordinateTransformY(X,Y,F,Ny)
%COORDINATETRANSFORM Summary of this function goes here
%   Detailed explanation goes here
% Maimouna Bocoum 13/02/2015

y_max = max(Y(:));
y_min = min(Y(:));

y = linspace(y_min,y_max,Ny);

%Fout = zeros(Ny,size(F,2));

 for i = 1:size(F,2)
 F_out = interp1(Y(:,i),F(:,i),y,'linear',0);
 Fout(:,i) = F_out;
 end

x = X(1,:);

[Xq,Yq] = meshgrid(x,y);

% figure;
% subplot(1,2,1)
% surf(X,Y,abs(F))
% shading interp
% view(0,90)
%   subplot(1,2,2)
%  % plot(Y(:,1),abs(F(:,1)))
%   surf(Xq,Yq,abs(Fout))
%   shading interp
% view(0,90)

end

