function []=Debug_interp_SPLINE()

%====================================================================
% Debugging the function interp_SPLINE
%
% Synopsis: []=Debug_interp_SPLINE
% 
% Written by Miki Elad on March 20th, 2005.
%====================================================================

cc=2;
if cc==1, % the sine function
    X=-pi:0.7:pi;
    Y=sin(X);   
    DY=cos(X);
    DDY=-sin(X);
    XX=-pi+0.7:0.001:pi-0.7;
else, % a polynomial
    X=-10:2:10;
    Y=0.0004*X.^5+0.2*X.^2+1;
    DY=0.002*X.^4+0.4*X;
    DDY=0.008*X.^3+0.4;
    XX=-9:0.01:9;
end;

YY1=interp1(X,Y,XX,'nearest'); % Matlab's interpolation
YY2=interp1(X,Y,XX,'linear'); % Matlab's interpolation
YY3=interp1(X,Y,XX,'spline'); % Matlab's interpolation
YY4=interp_SPLINE(X,Y,DY,DDY,XX);

figure(1); clf;
plot(XX,YY1,'r'); hold on;
plot(XX,YY2,'c'); 
plot(XX,YY3,'m'); 
plot(XX,YY4,'b');
plot(X,Y,'g.'); 
if cc==1,
    plot(XX,sin(XX),'y');
    disp(mean(abs(YY1-sin(XX)))); % Matlab's interpolation error
    disp(mean(abs(YY2-sin(XX)))); % Matlab's interpolation error
    disp(mean(abs(YY3-sin(XX)))); % Matlab's interpolation error
    disp(mean(abs(YY4-sin(XX)))); % our method error
else,
    YYt=0.0004*XX.^5+0.2*XX.^2+1;
    plot(XX,YYt,'y');
    disp(mean(abs(YY1-YYt))); % Matlab's interpolation error
    disp(mean(abs(YY2-YYt))); % Matlab's interpolation error
    disp(mean(abs(YY3-YYt))); % Matlab's interpolation error
    disp(mean(abs(YY4-YYt))); % our method error
end;    

return;
