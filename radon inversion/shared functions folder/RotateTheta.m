function [Iout,MMcorner] = RotateTheta(X,Z,Iin,theta,C)
% created by maimouna bocoum 
% 25/10/2017

% rotation with respect to box center in x , using fixed point of rotation
% C

X = X - C(1);
Z = Z - C(2);

% matrix of rotation (0,0) : center of rotation
M = [cos(theta), sin(theta) ; -sin(theta), cos(theta)];
% M = [cos(theta), -sin(theta) ; sin(theta), cos(theta)];

% rotation operate on all points
MM = M*[X(:)';Z(:)'] ;

% rotation of top-left corner
Mcorner =  M*[ min(X(:)) ; min(Z(:))  ] ;
MMcorner = Mcorner - [ min(X(:)) ; min(Z(:))  ] ;

% extract X,Z rotated coordinate
Xout = MM(1,:);
Zout = MM(2,:);

% reshape variable into original size
Xout = reshape(Xout,[size(X,1),size(X,2)]);
Zout = reshape(Zout,[size(X,1),size(X,2)]);

% map over previous coordinates
Iout = interp2(X,Z,Iin,Xout,Zout,'linear',0) ;

end

