% Jerome Faure% May 5 1997% This function calculates the derivative of a vector% V is the vector to derive% w is the variable which respect the derivation is made% h is the sample intervalfunction dd = derive(y,x)y = y(:)';x = x(:)';% interpolation dans un vecteur bien mieux r�solu :xmin = min(x);xmax = max(x);X = linspace(xmin,xmax,1e5);Y = interp1(x,y,X,'spline','extrap');h = X(2)-X(1);D = diff(Y)/h;DD = [D(1) D]';dd = interp1(X,DD,x,'nearest','extrap');