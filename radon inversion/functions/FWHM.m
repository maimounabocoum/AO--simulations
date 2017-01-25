function fwhm = FWHM(I,x)

% not : Et must be real
% Et must be the same length as t

X = linspace(min(x),max(x),200*length(x));

II = interp1(x,I,X,'linear','extrap');



MaxX = max(II);
X = X( II > MaxX/2);
fwhm = max(X) - min(X);

% figure()
% plot(x,I)
% hold on 
% plot(X,II(II > MaxX/2),'or')


