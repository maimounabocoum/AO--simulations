function []=Debug_XPolar_FFT

%======================================================================
% This function checks the errors in the polar FFT by testing the errors
% for the pseudo-polar, the S-pseudo-polar, and finally the X-pseudo-polar.
% 
% Synopsis: Debug_XPolar_FFT
% 
% Written by Michael Elad on March 20th, 2005.
%======================================================================

% Step 1: Testing the pseudo-polar accuracy

N=16; 
[Xc,Yc]=Create_Grid('D',[N,pi],'');
Trans=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the exact Transform matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=Brute_Force_Transform(x,Xc,Yc);
        Trans(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Approx_PP=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the fast Polar approximation matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=RectoPolar_Trans_New(x);
        % y=XPolar_Transform(x,3,OS,1,OS);
        Approx_PP(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Err_mat=Trans-Approx_PP; Y1=zeros(4*N^2,1);
for k=1:1:size(Err_mat,1),
    vec=Err_mat(k,:);
    if norm(vec)<1e-10, x=vec*0; else, x=vec/norm(vec); end;
    Y1(k)=abs(vec*x');
end;
figure(1); clf; plot(Y1);

%====> Conclusion - the Pseudo-Polar is exact and cannot account for the
%           strange behavior of the XPolar

% Step 2: Testing the SPolar accuracy
N=8; 
[Xc,Yc]=Create_Grid('S',[N,pi],'');
Trans=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the exact Transform matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=Brute_Force_Transform(x,Xc,Yc);
        Trans(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Approx_SP=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the fast Polar approximation matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=SPolar_Transform(x,4,4);
        Approx_SP(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Err_mat=Trans-Approx_SP; Y1=zeros(4*N^2,1);
for k=1:1:size(Err_mat,1),
    vec=Err_mat(k,:);
    if norm(vec)<1e-10, x=vec*0; else, x=vec/norm(vec); end;
    Y1(k)=abs(vec*x');
end;
figure(1); clf; plot(Y1);

delta=0.05;
data=[Xc(:),Yc(:),Y1];
L=length(-pi:delta:pi);
Result1=zeros(L,L);
ch=1; 
h=waitbar(0,'sweeping through 2D locations');
for hh=-pi:delta:pi,
    waitbar(ch/L);
    cv=1;
    for vv=-pi:delta:pi, 
        dist=(hh-data(:,1)).^2+(vv-data(:,2)).^2;
        pos=find(dist==min(dist));
        pos=pos(1);
        Result1(cv,ch)=data(pos,3);
        cv=cv+1;
    end;
    ch=ch+1;
end;
close(h);
figure(1); clf; 
imagesc(Result1); colorbar; axis off; axis image; 

%====> Conclusion - the S-Pseudo-Polar is behaving well, showing almost radial error
%           behavior as expected. Thus, we suspect the change from Spolar to
%           XPolar to cause the strange error structure

% Step 3: Testing the XPolar accuracy
N=8; 
[Xc,Yc]=Create_Grid('X',[N,pi],'');
Trans=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the exact Transform matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=Brute_Force_Transform(x,Xc,Yc);
        Trans(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Approx_XP=zeros(N^2*4,N^2);
count=1;
h=waitbar(0,'Building the fast Polar approximation matrix');
for k=1:1:N,
    waitbar(k/N);
    for j=1:1:N,
        x=zeros(N,N);
        x(k,j)=1;
        y=XPolar_Transform(x,3,4,1,20);
        Approx_XP(:,count)=y(:);
        count=count+1;
    end;
end;
close(h);

Err_mat=Trans-Approx_XP; Y1=zeros(4*N^2,1);
for k=1:1:size(Err_mat,1),
    vec=Err_mat(k,:);
    if norm(vec)<1e-10, x=vec*0; else, x=vec/norm(vec); end;
    Y1(k)=abs(vec*x');
end;
figure(1); clf; plot(Y1);

delta=0.05;
data=[Xc(:),Yc(:),Y1];
L=length(-pi:delta:pi);
Result1=zeros(L,L);
ch=1; 
h=waitbar(0,'sweeping through 2D locations');
for hh=-pi:delta:pi,
    waitbar(ch/L);
    cv=1;
    for vv=-pi:delta:pi, 
        dist=(hh-data(:,1)).^2+(vv-data(:,2)).^2;
        pos=find(dist==min(dist));
        pos=pos(1);
        Result1(cv,ch)=data(pos,3);
        cv=cv+1;
    end;
    ch=ch+1;
end;
close(h);
figure(2); clf; 
imagesc(Result1); colorbar; axis off; axis image; 

%====> Conclusion - the X-Pseudo-Polar seems to work well when the
%           OS2>>OS1.


return;
