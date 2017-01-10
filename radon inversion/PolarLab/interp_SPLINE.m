function YY=interp_SPLINE(X,Y,DY,DDY,XX)

%=====================================================================
% 1-D interpolation using function values and its 1st and 2nd derivatives. This function
% interpolates to find YY, the values of the underlying function Y at the points in the vector
% XX. The vector X specifies the points at which the data Y is given. A
% polynomial of order 6 is defined between each knots, and its 6 parameters
% are chosen based on the constraints from the data.
%
% Synopsis: YY=interp_SPLINE(X,Y,DY,DDY,XX)
%
% Inputs:  X      - absisa values where the value of the function are known
%               Y      - function values at X (i.e. f(x))
%              DX    - the function's derivative values at X (i.e. f'(x))
%               DDX  - the function sec ond derivatives as X (i.e. f''(x))
%               XX    - absisa locations where the function is desired.
%
% Outputs: YY     - interpolated valued at locations XX
%
% For debuging: replace all this code with YY=interp1(X,Y,XX,'spline');
%
% Written by Michael Elad on March 20th, 2005.
%=====================================================================

if X(2)-X(1)<0,
    X=fliplr(X);
    Y=fliplr(Y);
    DY=fliplr(DY);
    DDY=fliplr(DDY);
end;
Y(end+1)=0;
DY(end+1)=0;
DDY(end+1)=0;

method=1;
if method==1,

    Mat=[   1         0         0         0         0         0;
                 0         0         1         0         0        0;
                 0         0         0         0      0.5        0;
              -10       10        -6        -4     -1.5      0.5;
              15      -15         8         7       1.5       -1;
              -6          6       -3        -3      -0.5      0.5];

    Delta=X(2)-X(1);
    DY=DY*Delta;
    DDY=DDY*Delta^2;
    N=length(XX);
    YY=XX;
    % Pos=zeros(N,1);
    % for k=1:1:N, % sweeping through the desired values
    %     Pos(k)=max(find(XX(k)-X>=0));
    % end;

    % counting on X being uniformly sampled and including the origin, Pos can be computed
    % in an easier way. However, there is a dangerous inaccuracy involved
    Pos=floor((XX-X(1))/Delta+1e-8)+1;
    U=(XX-X(Pos))/Delta;
    Vec=[Y(Pos); Y(Pos+1); DY(Pos); DY(Pos+1); DDY(Pos); DDY(Pos+1)];
    YY=sum([U*0+1; U; U.^2; U.^3; U.^4; U.^5].*(Mat*Vec));

else, % using only first derivatives and not second
    Mat=   [ 1     0     0     0;
                 0     0     1     0;
                -3     3    -2    -1;
                 2    -2     1     1];

    Delta=X(2)-X(1);
    DY=DY*Delta;

    N=length(XX);
    YY=XX;
    Pos=floor((XX-X(1))/Delta+1e-8)+1;

    U=(XX-X(Pos))/Delta;
    Vec=[Y(Pos); Y(Pos+1); DY(Pos); DY(Pos+1)];
    YY=sum([U*0+1; U; U.^2; U.^3;].*(Mat*Vec));

end;

return;

%=====================================================
%                        The original code before the speeding changes
%=====================================================
if X(2)-X(1)<0,
    X=fliplr(X);
    Y=fliplr(Y);
    DY=fliplr(DY);
    DDY=fliplr(DDY);
end;

Mat=[  1         0         0         0         0         0;
    0         0         1         0         0        0;
    0         0         0         0      0.5        0;
    -10       10        -6        -4     -1.5      0.5;
    15      -15         8         7       1.5       -1;
    -6          6       -3        -3      -0.5      0.5];
% The inversion of the matrix:
% [1 0 0 0 0 0; 1 1 1 1 1 1; 0 1 0 0 0 0; 0 1 2 3 4 5; 0 0 2 0 0 0; 0 0 2 6 12 20]

Delta=X(2)-X(1);
DY=DY*Delta;
DDY=DDY*Delta^2;

N=length(XX);
YY=XX;
for k=1:1:N, % sweeping through the desired values
    Pos=max(find(XX(k)-X>=0));
    if XX(k)==X(Pos),
        YY(k)=Y(Pos);
    else,
        u=(XX(k)-X(Pos))/Delta;
        Vec=[Y(Pos); Y(Pos+1); DY(Pos); DY(Pos+1); DDY(Pos); DDY(Pos+1)];
        YY(k)=[1 u u^2 u^3 u^4 u^5]*Mat*Vec;
    end;
end;

return;
