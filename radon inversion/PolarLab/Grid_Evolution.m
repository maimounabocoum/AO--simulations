function [GridX1,GridY1,GridX2,GridY2,GridX3,GridY3]=Grid_Evolution(N,mode)

%=====================================================================
% This function describes the changes we perform when moving from the pseudo-polar 
% grid to the equi-spaced angles pseudo-polar grid, and from there to the Polar grid.
% 
% Synopsis: [GridX1,GridY1,GridX2,GridY2,GridX3,GridY3]=Grid_Evolution(N,mode)
% 
% Inputs:    N - Number of elements on each quadrant curve. The final grids will have 2N 
%                        rays with 2N elements on each.
%                  mode - 1 shows the grid ray by ray, and 0 shows it in one step
%
% Outputs:   GridX1,GridY1 - the pseudopolar grid
%                   GridX2,GridY2 - uniform angles instead of uniform distances along lines
%                   GridX3,GridY3 - uniform angles & uniform distances along rays
%
% Example:  [X1,Y1,X2,Y2,X3,Y3]=Grid_Evolution(16,1);
%
% Written by Michael Elad on Mrch 20th, 2005.
%=====================================================================


%-------------------------------------------------------------------------------------------------------
%                                 Stage 1 - Creating the pseudo-Polar grid
%-------------------------------------------------------------------------------------------------------

VX=zeros(2*N,N);
VY=zeros(2*N,N);
ll=-N:1:N-1;
theta_y=pi*ll/N;
m=-N/2:1:N/2-1;
for k=1:1:length(theta_y),
    theta_x=2*m*theta_y(k)/N;    
    VX(k,:)=theta_x;
    VY(k,:)=theta_y(k)*ones(1,N);
end;

HX=zeros(2*N,N);
HY=zeros(2*N,N);
ll=-N:1:N-1;
theta_x=pi*ll/N;
m=-N/2+1:1:N/2;
for k=1:1:2*N,
    theta_y=2*m*theta_x(k)/N;    
    HX(k,:)=theta_x(k)*ones(1,N);
    HY(k,:)=theta_y;
end;

GridX1=[fliplr(VX),HX]; % we have to flip in order to get a smoothly rotating rays
GridY1=[fliplr(VY),HY];

figure(1); clf; 
if mode==1,
    plot(0,0,'w.'); 
    xlabel('The starting grid: Pseudo (recto)-Polar');
    hold on;
    axis equal;
    axis([-pi pi -pi pi]);
    for k=1:1:2*N,
        plot(GridX1(:,k),GridY1(:,k),'b.','Markersize',10);
        pause;
    end;
else,
    plot(GridX1,GridY1,'b.','Markersize',10);    
    axis equal;
    axis([-3.5 3.5 -3.5 3.5]);
end;

%-------------------------------------------------------------------------------------------------------
%                                 Stage 2 - Modifying to obtain uniform angles
%-------------------------------------------------------------------------------------------------------

VX=zeros(2*N,N);
VY=zeros(2*N,N);
ll=-N:1:N-1;
theta_y=pi*ll/N;
m=-N/2:1:N/2-1;
for k=1:1:length(theta_y),
    theta_x=theta_y(k)*tan(m*pi/N/2);    
    VX(k,:)=theta_x;
    VY(k,:)=theta_y(k)*ones(1,N);
end;

HX=zeros(2*N,N);
HY=zeros(2*N,N);
ll=-N:1:N-1;
theta_x=pi*ll/N;
m=-N/2+1:1:N/2;
for k=1:1:2*N,
    theta_y=theta_x(k)*tan(m*pi/N/2);    
    HX(k,:)=theta_x(k)*ones(1,N);
    HY(k,:)=theta_y;
end;

GridX2=[fliplr(VX),HX]; % we have to flip in order to get a smoothly rotating rays
GridY2=[fliplr(VY),HY];

figure(2); clf; 
if mode==1,
    plot(0,0,'w.'); 
    xlabel('The second grid: Rays are rotated to give equiangular spread');
    hold on;
    axis equal;
    axis([-pi pi -pi pi]);
    for k=1:1:2*N,
        plot(GridX2(:,k),GridY2(:,k),'b.','Markersize',10);
        pause;
    end;
else,
    plot(GridX2,GridY2,'b.','Markersize',10);    
    axis equal;
    axis([-3.5 3.5 -3.5 3.5]);
end;

%-------------------------------------------------------------------------------------------------------
%                        Stage 3 - Modifying to obtain uniform distances on all rays
%-------------------------------------------------------------------------------------------------------

GridX3=GridX2;
GridY3=GridY2;

for k=1:1:2*N,
    Factor=sqrt(GridX3(k,:).^2+GridY3(k,:).^2);
    if min(Factor>0), 
        Factor=Factor./min(Factor);
        GridX3(k,:)=GridX3(k,:)./Factor;
        GridY3(k,:)=GridY3(k,:)./Factor;    
    end;
end;

figure(3); clf; 
if mode==1,
    plot(0,0,'w.'); 
    xlabel('The final grid: Rays are shrinked to become pure polar grid');
    hold on;
    axis equal;
    axis([-pi pi -pi pi]);
    for k=1:1:2*N,
        plot(GridX3(:,k),GridY3(:,k),'b.','Markersize',10);
        pause;
    end;
else,
    plot(GridX3,GridY3,'b.','Markersize',10);    
    axis equal;
    axis([-3.5 3.5 -3.5 3.5]);
end;

return;   