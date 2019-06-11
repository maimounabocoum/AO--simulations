function [x_out y_out I_out] = AffineTransform(x,y,I,A,B)
%AFFINETRANSFORM Summary of this function goes here
% (x',y') = A(x,y)+ B in the matrix relation
% we know f(x,y) = I
% and we return g(x,y) = f(A*(x,y)+B)
%   Detailed explanation goes here
x_out = x;
y_out = y;

[X Y] = meshgrid(x,y);
xx = 1/det(A)*( A(2,2)*(X-B(1)) -A(1,2)*(Y -B(2)));

yy = 1/det(A)*(-A(2,1)*(X-B(1))+ A(1,1)*(Y -B(2)));

I_out = interp2(X,Y,I,xx,yy,'linear',0);

end

